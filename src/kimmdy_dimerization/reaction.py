from kimmdy.recipe import (
    Bind,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    get_atomnrs_from_plumedid,
)
import logging
from kimmdy.parsing import (
    read_plumed,
    read_distances_dat,
)

import MDAnalysis as mda

logger = logging.getLogger("kimmdy.dimerization")


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
        length_of_traj = len(universe.trajectory)
        end_time = universe.trajectory[-1].time
        logger.info(f"Length of trajectory is {length_of_traj}, final time is {end_time} ps.")
        logger.debug("Getting recipe for reaction: Dimerization")

        recipes = []
        # check distances file for values lower than the cutoff

        recipes.append(
                    Recipe(
                        recipe_steps=[
                            Bind(atom_id_1="1", atom_id_2="2"),
                        ],
                        rates=[1],
                        timespans=[(0, end_time)],
                    )
                )

        return RecipeCollection(recipes)