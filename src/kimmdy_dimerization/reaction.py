from kimmdy.recipe import (
    Bind,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
import logging

from MDAnalysis.analysis.dihedrals import Dihedral
import MDAnalysis as mda
import numpy as np
import math
from itertools import combinations

logger = logging.getLogger("kimmdy.dimerization")


def calculate_rate(k1_in, k2_in, d0_in, n0_in, distance_in, angle_in):
    return math.exp(-(k1_in * abs(distance_in - d0_in) + k2_in * abs(angle_in - n0_in)))


class DimerizationReaction(ReactionPlugin):
    """A Reaction Plugin for Dimerization in DNA
    """

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

        recipes = []
        time_start = 0
        for frame_idx, rates in enumerate(rates_time_resolved):
            if frame_idx != len(universe.trajectory) - 1:
                time_end = universe.trajectory[frame_idx + 1].time
            else:
                time_end = time_start
            if time_start == 0:
                for rate in rates:
                    recipes.append(
                        Recipe(
                            recipe_steps=[
                                Bind(atom_id_1="15", atom_id_2="46"),
                                Bind(atom_id_1="12", atom_id_2="44")
                            ],
                            rates=[rate[2]],
                            timespans=[(time_start, time_end)],
                        )
                    )
            else:
                for rate in rates:
                    recipes.append(
                        Recipe(
                            recipe_steps=[
                                Bind(atom_id_1="14", atom_id_2="46"),
                                Bind(atom_id_1="12", atom_id_2="44")
                            ],
                            rates=[rate[2]],
                            timespans=[(time_start, time_end)],
                        )
                    )
            time_start = time_end

        logger.info(f"Found rates {rates_time_resolved}")
        return RecipeCollection(recipes)
