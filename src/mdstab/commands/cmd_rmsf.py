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
"""Subcommand to calculate the r.m.s.f."""
from pathlib import Path

import click_extra as click
import MDAnalysis as mda
import numpy as np
from click_extra.logging import logger as clog
from MDAnalysis.analysis import rms

from .. import _MASK
from .. import PathLike
from .. import __copyright__
from .. import create_logger


_help = f"{__copyright__}\n\nCalculate the r.m.s.f. of a trajectory."
_short_help = "Calculate r.m.s.f."


@click.command("rmsf", short_help="Calculate root mean square fluctuations.")
@click.option(
    "-s",
    "--topology",
    metavar="FILE",
    default="amber.prmtop",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    help="Topology file",
)
@click.option(
    "-f",
    "--trajectory",
    metavar="FILE",
    default="amber.pdb",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    help="Trajectory file",
)
@click.option(
    "-o",
    "--outfile",
    metavar="FILE",
    default=Path.cwd() / "rmsf.npy",
    type=click.Path(exists=False, dir_okay=False, resolve_path=True, path_type=Path),
    help="r.m.s.f. data file",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    default=Path.cwd() / "rmsf.log",
    type=click.Path(exists=False, dir_okay=False, resolve_path=True, path_type=Path),
    help="Log file",
)
@click.option(
    "-b",
    "start",
    metavar="START",
    default=1,
    type=click.IntRange(min=1, clamp=True),
    help="Starting trajectory frame",
)
@click.option(
    "-e",
    "stop",
    metavar="STOP",
    default=0,
    type=click.IntRange(min=0, clamp=True),
    help="Final trajectory frame",
)
@click.option(
    "--dt",
    "offset",
    metavar="OFFSET",
    default=1,
    type=click.IntRange(min=1, clamp=True),
    help="Trajectory output offset (0 = last frame)",
)
@click.option(
    "-t",
    "--type",
    "calc_type",
    metavar="TYPE",
    default="ca",
    type=click.Choice(_MASK.keys(), case_sensitive=True),
    help=f"Atom selection {list(_MASK.keys())}",
)
def cli(
    topology: PathLike,
    trajectory: PathLike,
    outfile: PathLike,
    logfile: PathLike,
    start: int,
    stop: int,
    offset: int,
    calc_type: str,
) -> None:
    """Calculate the root mean square fluctuations of both heavy and selected atoms."""
    """Prepare various Amber input files to run simulations."""
    click.echo(__copyright__)

    # Setup logging
    logger: clog = create_logger(__name__, clog.wrapped_logger.level, logfile=logfile)

    logger.debug(f"Loading universe using {topology} and {trajectory}")
    universe = mda.Universe(topology, trajectory)
    logger.debug(f"Selecting atoms with keywords '{_MASK[calc_type]}'")
    atoms: mda.AtomGroup = universe.select_atoms(_MASK[calc_type])
    end = universe.trajectory.n_frames if stop < 1 else stop

    logger.info(f"Calculating the r.m.s.f. for {atoms.n_atoms} atoms.")
    logger.debug(f"Start frame: {start}; stop frame: {end}; frame step: {offset}")
    rmsf = rms.RMSF(atoms).run(start=start - 1, stop=end, step=offset)

    logger.info(f"Saving r.m.s.f. data to {outfile}")
    np.save(outfile, rmsf.results.rmsf)
