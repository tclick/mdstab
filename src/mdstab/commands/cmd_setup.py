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
"""Create simulation and analysis subdirectories."""
from itertools import product
from pathlib import Path

import click_extra as click
from click_extra.logging import logger as clog

from .. import PathLike
from .. import __copyright__
from .. import create_logger


_help = f"{__copyright__}\n\nCreate subdirectories for molecular dynamic simulations."
_short_help = "Setup simulation and analysis subdirectories"


@click.command("setup", help=_help, short_help=_short_help)
@click.option(
    "-o",
    "--outdir",
    metavar="DIR",
    default=Path.cwd(),
    type=click.Path(file_okay=False, resolve_path=True, path_type=Path),
    help="Parent directory",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    default=Path.cwd() / "setup.log",
    type=click.Path(dir_okay=False, resolve_path=True, path_type=Path),
    help="Log file",
)
def cli(outdir: PathLike, logfile: PathLike) -> None:
    """Create subdirectories for molecular dynamic simulations."""
    click.echo(__copyright__)

    # Setup logging
    logger = create_logger(__name__, clog.wrapped_logger.level, logfile=logfile)

    dirs = ("Prep", "Equil", "Prod", "Analysis")
    for _ in dirs:
        directory = Path(outdir) / _
        logger.info(f"Creating {directory}")
        directory.mkdir(mode=0o755, parents=True, exist_ok=True)

    equil = ("min", "md")
    subsection = (1, 2, 11, 12, 13, 14, 15, 16)
    for x, y in product(equil, subsection):
        directory = Path(outdir) / "Equil" / f"{x}{y:d}"
        if directory != "min16":
            logger.info(f"Creating {directory}")
            directory.mkdir(mode=0o755, parents=True, exist_ok=True)

    prod = ("initial", "production")
    for _ in prod:
        directory = Path(outdir) / "Prod" / _
        logger.info(f"Creating {directory}")
        directory.mkdir(mode=0o755, parents=True, exist_ok=True)
