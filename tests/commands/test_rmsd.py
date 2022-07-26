#  mdstab — Molecular Dynamics Setup and Trajectory Analysis for Biomolecules
#  Copyright (C) 2022  Timothy H. Click, Ph.D.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Test the r.m.s.d. subcommand."""
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from mdstab import cli

from ..datafile import TOPWW
from ..datafile import TRJWW


class TestRmsd:
    """Test rmsd subcommand."""

    @pytest.fixture
    def runner(self) -> CliRunner:
        """Fixture for invoking command-line interfaces.

        Returns
        -------
        CliRunner
            Test runner
        """
        return CliRunner()

    def test_help(self, runner: CliRunner) -> None:
        """Test help output.

        Parameters
        ----------
        runner : CliRunner
            Test runner
        """
        result = runner.invoke(cli.main, ["rmsd", "-h"])
        assert result.exit_code == 0
        assert "--verbosity" in result.output

    def test_rmsd(self, runner: CliRunner, tmp_path: Path) -> None:
        """Test rmsd subcommand.

        Parameters
        ----------
        runner : CliRunner
            Test runner
        tmp_path : Path
            temporary directory
        """
        logfile = tmp_path / "rmsd.log"
        outfile = tmp_path / "rmsd.npy"
        result = runner.invoke(
            cli.main,
            [
                "rmsd",
                "-s",
                TOPWW,
                "-f",
                TRJWW,
                "-o",
                outfile.as_posix(),
                "-l",
                logfile.as_posix(),
            ],
        )

        assert result.exit_code == 0
        assert logfile.is_file()
        assert outfile.exists()

        rmsd = np.load(outfile)
        assert np.all(rmsd[:, 2] > 0)
