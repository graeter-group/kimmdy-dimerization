"""File for pytest configuration in python and fixture definition.

The name 'conftest.py' is recognized by pytest to execute it before tests.
"""

import pytest
import shutil
import os
from pathlib import Path


## fixtures for setup and teardown ##
@pytest.fixture
def arranged_tmp_path(tmp_path: Path, request: pytest.FixtureRequest):
    """Arrange temporary directory for tests.

    With files for the test and a symlink to forcefield.
    """

    # if fixture was parameterized, use this for directory with input files
    if hasattr(request, "param"):
        file_dir = Path(__file__).parent / request.param
    # else use stem of requesting file to find directory with input files
    else:
        file_dir = Path(__file__).parent / request.path.stem
    # arrange tmp_path

    shutil.copytree(file_dir, tmp_path, dirs_exist_ok=True)
    assetsdir = Path(__file__).parent / "assets"

    if not (tmp_path / "amber14sb_OL21_mod.ff").exists():
        Path(tmp_path / "amber14sb_OL21_mod.ff").symlink_to(
            assetsdir / "amber14sb_OL21_mod.ff",
            target_is_directory=True,
        )
    # change cwd to tmp_path
    os.chdir(tmp_path.resolve())
    return tmp_path
