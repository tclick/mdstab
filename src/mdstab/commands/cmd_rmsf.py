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
import seaborn as sns
from click_extra.logging import logger as clog
from matplotlib import pyplot as plt
from MDAnalysis.analysis import rms
from modin import pandas as pd

from .. import _MASK
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
@click.option(
    "--label",
    metavar="LABEL",
    default=10,
    show_default=True,
    type=click.IntRange(min=1, clamp=True),
    help="Spacing for tick labels",
)
@click.option("--image", is_flag=True, help="Save correlation heatmap")
def cli(
    topology: Path,
    trajectory: Path,
    outfile: Path,
    logfile: Path,
    start: int,
    stop: int,
    offset: int,
    calc_type: str,
    label: int,
    image: bool,
) -> None:
    """Calculate the root mean square fluctuations of both heavy and selected atoms."""
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

    if image:
        filename = outfile.with_suffix(".png")

        fig: plt.Figure = plt.figure(figsize=plt.figaspect(1.0))
        title = _MASK[calc_type]
        df = pd.DataFrame(
            [atoms.resnums, rmsf.results.rmsf], columns=["Residue", "r.m.s.f. (Å)"]
        )
        df[df.columns[0]] = df[df.columns[0]].astype(int)
        ax = fig.add_subplot(1, 1, 1)
        sns.barplot(x="Residue", y="r.m.s.f. (Å)", data=df, color="blue", ax=ax)
        ax.set_title(f"{title[i]}")
        ax.set_xticks(np.arange(label - 1, rmsf.results.rmsf.size, label))
        fig.suptitle("Root mean square fluctuations")
        logger.info(f"Saving image to {filename}")
        fig.autofmt_xdate(rotation=45)
        fig.tight_layout()
        fig.savefig(filename, dpi=600)
