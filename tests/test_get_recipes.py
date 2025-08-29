from pathlib import Path
import numpy as np
from kimmdy.recipe import Bind, CustomTopMod, Relax
import pytest
from kimmdy_dimerization.reaction import DimerizationReaction


def test_get_recipe_collection(mocker, arranged_tmp_path):
    m_runmanager = mocker.Mock()
    m_config = mocker.Mock()
    m_config.k1 = 2.017017017017017
    m_config.k2 = 0.03003003003003003
    m_config.d0 = 0.157177
    m_config.n0 = 16.743651884789273
    m_config.reslist = "all"

    m_files = mocker.Mock()
    m_files.input = {
        "gro": Path("equilibrium_short.gro"),
        "trr": Path("equilibrium_short.trr"),
    }
    m_files.outputdir = arranged_tmp_path

    m_runmanager.config.reactions.test = m_config
    dimerization = DimerizationReaction("test", m_runmanager)

    recipe_collection = dimerization.get_recipe_collection(m_files)
    assert (arranged_tmp_path / "reaction_rates.csv").exists()
    assert len(recipe_collection.recipes) == 1

    # Test all recipe steps
    assert recipe_collection.recipes[0].recipe_steps[0] == Bind(
        atom_ix_1=13, atom_ix_2=45
    )
    assert recipe_collection.recipes[0].recipe_steps[1], Bind(
        atom_ix_1=11, atom_ix_2=43
    )

    def _change_top():
        pass

    assert (
        recipe_collection.recipes[0]
        .recipe_steps[2]
        .__almost_eq__(CustomTopMod(_change_top))
    )
    assert recipe_collection.recipes[0].recipe_steps[3] == (Relax())

    ts = np.array(recipe_collection.recipes[0].timespans)
    a = np.arange(11)
    b = np.arange(11) + 1
    b[-1] = 10
    np.testing.assert_array_equal(np.stack((a, b), 1), ts)

    assert len(recipe_collection.recipes[0].rates) == 11
    assert pytest.approx(min(recipe_collection.recipes[0].rates), 0.01) == 0.219
    assert pytest.approx(max(recipe_collection.recipes[0].rates), 0.01) == 0.559
