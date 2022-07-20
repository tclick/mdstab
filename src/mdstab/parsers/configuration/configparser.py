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
"""Parse a configuration file."""
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import List

import click
import pylibyaml  # noqa: F401
import yaml

from ... import PathLike


logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class Config:
    """Configuration data."""

    def update(self, **kwargs: Any) -> None:
        """Add or update attributes.

        Parameters
        ----------
        kwargs : Any
            Dict of attributes to add
        """
        for key, value in kwargs.items():
            setattr(self, key, value)


def configure(ctx: click.Context, param: List[Any], filename: PathLike) -> None:
    """Read a configuration file and pouplate click.Context.default_map with contents.

    Parameters
    ----------
    ctx : click.Context
        Context object
    param : list
        List of parameters
    filename : PathLike
        configuration file
    """
    if Path(filename).exists():
        logger.debug(f"Using configuration file: {filename}")
        with open(filename) as config_file:
            ctx.default_map = yaml.safe_load(config_file)
