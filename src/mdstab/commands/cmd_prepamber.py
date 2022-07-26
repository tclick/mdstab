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
"""Prepare Amber simulation files."""
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import click_extra as click
import MDAnalysis as mda
from click_extra.logging import logger as clog
from jinja2 import Environment
from jinja2 import PackageLoader

from .. import PathLike
from .. import __copyright__
from .. import create_logger


_help = f"{__copyright__}\n\nPrepare Amber simulation files."
_short_help = "Prepare Amber input files."


@dataclass
class Data:
    """Data for jinja2 template."""

    temp1: float
    temp2: float
    res0: int
    res1: int
    ions0: int
    ions1: int
    solvent0: int
    solvent1: int
    force: float
    simdir: PathLike
    prefix: str
    amberhome: PathLike
    pmemd: str


def write_template(
    env: Environment, temploc: str, data: Data, simdir: PathLike, logger: clog
) -> None:
    """Write jinja2 template to a file.

    Parameters
    ----------
    env : Environment
        jinja2 environment
    temploc : str
        template location
    data : Data
        template data
    simdir : PathLike
        subdirectory for files
    logger : clog
        logger
    """
    template_dir = Path("templates") / "amber" / f"{temploc}"
    env.loader = PackageLoader("mdstab", package_path=template_dir.as_posix())
    for template_file in env.loader.list_templates():
        subdir = (
            Path(simdir) / temploc.title() / Path(template_file).stem
            if temploc != "scripts"
            else Path(simdir) / temploc.title()
        )
        subdir.mkdir(parents=True, exist_ok=True)
        filename = subdir / template_file
        input_file = (
            filename.with_suffix(".sh")
            if temploc == "scripts"
            else filename.with_suffix(".in")
        )
        with open(input_file, mode="w") as inf:
            template = env.get_template(template_file)
            logger.info(f"Writing script to {input_file}")
            print(template.render(data=data), file=inf)


@click.command("prepamber", help=_help, short_help=_short_help)
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
    "-d",
    "--simdir",
    metavar="DIR",
    default=Path.cwd(),
    type=click.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Simulation subdirectory",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PREFIX",
    default=Path.cwd().stem,
    help="Prefix for various output files",
)
@click.option(
    "--temp1",
    metavar="TEMP",
    default=100.0,
    type=click.IntRange(min=1.0, clamp=True),
    help="Initial temperature (K)",
)
@click.option(
    "--temp2",
    metavar="TEMP",
    default=300.0,
    type=click.IntRange(min=1.0, clamp=True),
    help="Final temperature (K)",
)
@click.option(
    "--force",
    metavar="FORCE",
    default=100.0,
    type=click.IntRange(min=1.0, clamp=True),
    help="Restraint force (kcal/mol/A^2",
)
@click.option(
    "-l",
    "--logfile",
    metavar="LOG",
    default=Path.cwd() / "prepare.log",
    type=click.Path(exists=False, resolve_path=True, path_type=Path),
    help="Log file",
)
@click.option(
    "--amberhome",
    metavar="DIR",
    type=click.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Location of Amber files",
)
@click.option(
    "--type",
    "outtype",
    default="all",
    type=click.Choice("equil prod shell all".split(), case_sensitive=False),
    help="Which output files to create",
)
def cli(
    topology: Path,
    trajectory: Path,
    simdir: Path,
    prefix: str,
    temp1: float,
    temp2: float,
    force: float,
    logfile: Path,
    amberhome: Path,
    outtype: str,
) -> None:
    """Prepare various Amber input files to run simulations."""
    click.echo(__copyright__)

    # Setup logging
    logger = create_logger(__name__, clog.wrapped_logger.level, logfile=logfile)

    universe = mda.Universe(topology, trajectory)
    protein = universe.select_atoms("protein")
    ions = universe.select_atoms("name Na+ K+ Cl-")
    solvent = universe.select_atoms("resname WAT")
    try:
        amberhome = (
            Path(os.environ["AMBERHOME"]) if amberhome is None else Path(amberhome)
        )
    except KeyError:
        logger.exception("AMBERHOME environment variable not defined.", exc_info=True)

    simdir = Path(simdir)
    data = Data(
        temp1=temp1,
        temp2=temp2,
        res0=protein.resnums[0],
        res1=protein.resnums[-1],
        ions0=ions.resnums[0],
        ions1=ions.resnums[-1],
        solvent0=solvent.resnums[0],
        solvent1=solvent.resnums[-1],
        force=force,
        simdir=Path(simdir),
        prefix=prefix,
        amberhome=amberhome,
        pmemd=(
            "pmemd.MPI"
            if shutil.which(amberhome / "bin" / "pmemd.MPI") is not None
            else "pmemd"
        ),
    )

    env = Environment(autoescape=True)
    temploc = []
    if outtype.lower() == "equil":
        temploc.append("equil")
    elif outtype.lower() == "prod":
        temploc.append("prod")
    elif outtype.lower() == "shell":
        temploc.append("scripts")
    else:
        temploc.extend(["equil", "prod", "scripts"])

    for _ in temploc:
        write_template(env, _, data, simdir, logger)

    # Write the shell scripts
    if "scripts" in temploc:
        scripts = Path(simdir) / "Scripts"
        for _ in scripts.glob("*.sh"):
            _.chmod(0o755)
