from kimmdy.recipe import (
    Bind,
    Recipe,
    RecipeCollection,
    CustomTopMod,
    Relax
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
import logging
from kimmdy.topology.topology import Topology
from MDAnalysis.analysis.dihedrals import Dihedral
import MDAnalysis as mda
import numpy as np
import math
import json
from importlib.resources import files
from itertools import combinations

logger = logging.getLogger("kimmdy.dimerization")


def calculate_rate(k1_in, k2_in, d0_in, n0_in, distance_in, angle_in):
    return math.exp(-(k1_in * abs(distance_in - d0_in) + k2_in * abs(angle_in - n0_in)))


class DimerizationReaction(ReactionPlugin):
    """A Reaction Plugin for Dimerization in DNA
    """
    @staticmethod
    def change_top(res_a, res_b):
        change_dict = {"C6": "CT", "C5": "CT", "H6": "H1", "N1": "N"}

        path = files("kimmdy_dimerization.data") / "new_charges.json"
        with path.open("r", encoding="utf-8") as f:
            new_charges = json.load(f)
        res_a = str(res_a)
        res_b = str(res_b)

        def f(top: Topology) -> Topology:
            # Determine newly bonded atoms
            for atom in top.atoms.values():
                if atom.resnr == res_a and atom.atom == "C5":
                    c5_a = atom
                if atom.resnr == res_a and atom.atom == "C6":
                    c6_a = atom
                if atom.resnr == res_b and atom.atom == "C5":
                    c5_b = atom
                if atom.resnr == res_b and atom.atom == "C6":
                    c6_b = atom

            # Find improper dihedrals at C5 and C6 that need to be removed
            dihedrals_to_remove = []
            for dihedral_key in top.improper_dihedrals.keys():
                if (c5_a.nr in dihedral_key and c6_a.nr in dihedral_key) or (
                        c5_b.nr in dihedral_key and c6_b.nr in dihedral_key):
                    dihedrals_to_remove.append(dihedral_key)
            for dihedral_key in dihedrals_to_remove:
                logger.info(f"Removed improper dihedral {dihedral_key}")
                top.improper_dihedrals.pop(dihedral_key, None)

            # Change residue types
            for atom in top.atoms.values():
                if atom.resnr == res_a or atom.resnr == res_b:
                    atom.residue = atom.residue.replace("T", "D")
            # Change atomtypes
            for atom in top.atoms.values():
                if atom.resnr == res_a or atom.resnr == res_b:
                    if atom.atom in change_dict.keys():
                        atom.type = change_dict[atom.atom]
            # Change charges (already present in ff for new residue type)
            for atom in top.atoms.values():
                if atom.resnr == res_a or atom.resnr == res_b:
                    atom.charge = new_charges[atom.residue][atom.atom]

            # Remove faulty pairs
            top.pairs.pop((c6_a.nr, c6_b.nr), None)
            top.pairs.pop((c5_a.nr, c5_b.nr), None)
            top.pairs.pop((c6_a.nr, c5_a.nr), None)
            top.pairs.pop((c6_b.nr, c5_b.nr), None)

            return top

        return CustomTopMod(f)

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger

        # Get values from config
        k1 = self.config.k1  # Distance scaling [1/nm]
        k2 = self.config.k2  # Angle scaling [1/deg]
        d0 = self.config.d0  # Optimal distance [nm]
        n0 = self.config.n0  # Optimal angle [deg]

        gro = files.input["gro"]
        trr = files.input["trr"]
        universe = mda.Universe(str(gro), str(trr))
        c5s = universe.select_atoms("name C5 and resname DT5 DT DT3")
        c6s = universe.select_atoms("name C6 and resname DT5 DT DT3")
        c5_c6s = [(c5, c6) for c5, c6 in zip(c5s, c6s)]
        residue_dict_c5 = {int(c5.resid): c5.ix for c5 in c5s}
        residue_dict_c6 = {int(c6.resid): c6.ix for c6 in c6s}

        # Dihedrals
        dihedrals_time_resolved = [[] for _ in range(0, len(universe.trajectory))]
        for reactive_four in combinations(c5_c6s, r=2):
            dihedral_group = mda.AtomGroup(
                [reactive_four[0][0], reactive_four[0][1], reactive_four[1][1], reactive_four[1][0]])
            dih = Dihedral([dihedral_group])
            dih.run()
            for time_idx, ang in enumerate(dih.results.angles):
                dihedrals_time_resolved[time_idx].append(
                    (int(reactive_four[0][0].resid), int(reactive_four[1][0].resid), float(ang[0])))

        # Distances
        dists_time_resolved = []
        for _ in universe.trajectory:
            vecs_c5_1_c5_2 = [(c5_1.resid, c5_2.resid, c5_2.position - c5_1.position) for c5_1, c5_2 in
                              combinations(c5s, r=2)]
            vecs_c6_1_c6_2 = [(c6_1.resid, c6_2.resid, c6_2.position - c6_1.position) for c6_1, c6_2 in
                              combinations(c6s, r=2)]
            dists = [(int(vec_c5_1_c5_2[0]), int(vec_c5_1_c5_2[1]),
                      0.1 * float(np.linalg.norm(0.5 * (vec_c5_1_c5_2[2] + vecs_c6_1_c6_2[2]))))
                     for vec_c5_1_c5_2, vecs_c6_1_c6_2 in zip(vecs_c5_1_c5_2, vecs_c6_1_c6_2)]
            dists_time_resolved.append(dists)

        # Rate calculation
        rates_time_resolved = [[] for _ in range(0, len(universe.trajectory))]
        for time_idx, (distances, angles) in enumerate(zip(dists_time_resolved, dihedrals_time_resolved)):
            for distance, angle in zip(distances, angles):
                rates_time_resolved[time_idx].append(
                    (distance[0], distance[1], calculate_rate(k1, k2, d0, n0, distance[2], angle[2])))

        # Group by reaction
        reactions = {}
        time_start = 0
        for frame_idx, rates in enumerate(rates_time_resolved):
            if frame_idx != len(universe.trajectory) - 1:
                time_end = universe.trajectory[frame_idx + 1].time
            else:
                time_end = time_start
            for rate in rates:
                if (rate[0], rate[1]) not in reactions.keys():
                    reactions[(rate[0], rate[1])] = [(rate, time_start, time_end)]
                else:
                    reactions[(rate[0], rate[1])].append((rate, time_start, time_end))
            time_start = time_end
        logger.info(f"Found reactions with rates {reactions}")

        recipes = []
        for reaction_key in reactions.keys():
            reaction = reactions[reaction_key]
            rates = [a[0][2] for a in reaction]
            timespans = [(a[1], a[2]) for a in reaction]
            res_a = reaction_key[0]
            res_b = reaction_key[1]

            steps = [
                Bind(atom_id_1=str(residue_dict_c5[res_a]+1), atom_id_2=str(residue_dict_c5[res_b]+1)),
                Bind(atom_id_1=str(residue_dict_c6[res_a]+1), atom_id_2=str(residue_dict_c6[res_b]+1)),
                self.change_top(res_a, res_b),
                Relax()
            ]
            recipes.append(
                Recipe(
                    recipe_steps=steps,
                    rates=rates,
                    timespans=timespans
                )
            )
        logger.info(f"Recipe collection {recipes}")

        return RecipeCollection(recipes)
