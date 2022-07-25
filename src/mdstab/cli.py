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
"""Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mmdta` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``mdta.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``mdta.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import logging
import sys
from pathlib import Path
from typing import Any
from typing import List
from typing import Optional

import click_extra as click

from . import PathLike
from . import __copyright__
from . import __version__


CONTEXT_SETTINGS = dict(
    auto_envvar_prefix="COMPLEX", help_option_names=["-h", "--help"], show_default=True
)
logger = logging.getLogger()
cmd_folder = Path(__file__).parent.joinpath("commands").resolve()


class Environment:
    """Context manager for click command-line interface."""

    def log(self, msg: str, *args: List[str]) -> None:
        """Log a message to stderr.

        Parameters
        ----------
        msg : str
            Log message
        args : list of strings
            additional arguments to be passed
        """
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg: str, *args: List[str]) -> None:
        """Log a message to stderr only if verbose is enabled.

        Parameters
        ----------
        msg : str
            Log message
        args : list of strings
            additional arguments to be passed
        """
        if self.verbose:
            self.log(msg, *args)


pass_environment = click.make_pass_decorator(Environment, ensure=True)


class ComplexCLI(click.MultiCommand):
    """Complex command-line options with subcommands for fluctmatch."""

    def list_commands(self, ctx: click.Context) -> Optional[List[str]]:
        """List available commands.

        Parameters
        ----------
        ctx : `Context`
            click context

        Returns
        -------
            List of available commands
        """
        rv = []
        for filename in Path(cmd_folder).iterdir():
            if filename.name.endswith(".py") and filename.name.startswith("cmd_"):
                rv.append(filename.name[4:-3])
        rv.sort()
        return rv

    def get_command(self, ctx: click.Context, name: str) -> Optional[Any]:
        """Run the selected command.

        Parameters
        ----------
        ctx : `Context`
            click context
        name : str
            command name

        Returns
        -------
            The chosen command if present
        """
        try:
            mod = __import__(f"mdstab.commands.cmd_{name}", None, None, ["cli"])
        except ImportError:
            return None
        return mod.cli


@click.command(cls=ComplexCLI, context_settings=CONTEXT_SETTINGS, help=__copyright__)
@click.colorize.version_option(version=__version__)
@click.option(
    "--home",
    type=click.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
    help="Changes the folder to operate on.",
)
@click.config_option()
@pass_environment
def main(ctx: click.Context, home: PathLike) -> None:
    """Molecular dynamics setup and trajectory analysis for biomolecules."""
    if home is not None:
        ctx.home = home
