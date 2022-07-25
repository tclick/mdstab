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
"""Test setup subcommand."""
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from click.testing import CliRunner
from click_extra import platform

from mdstab import cli


class TestSetup:
    """Test setup subcommand."""

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
        result = runner.invoke(cli.main, ["setup", "-h"])
        assert result.exit_code == 0
        assert "--verbosity" in result.output

    def test_setup(self, runner: CliRunner) -> None:
        """Test setup of MD subdirectories.

        Parameters
        ----------
        runner : CliRunner
            Test runner
        """
        ignore_cleanup = platform.is_windows()
        with TemporaryDirectory(ignore_cleanup_errors=ignore_cleanup) as outdir:
            out_dir = Path(outdir)
            logfile = out_dir / "setup.log"
            result = runner.invoke(
                cli.main, ["setup", "-o", out_dir.as_posix(), "-l", logfile.as_posix()]
            )

            assert result.exit_code == 0
            assert logfile.is_file()
            assert (out_dir / "Prep").is_dir()
