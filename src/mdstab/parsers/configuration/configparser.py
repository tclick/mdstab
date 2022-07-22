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
from dataclasses import make_dataclass
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

import click
import pylibyaml  # noqa: F401
import yaml

from ... import PathLike


logger: logging.Logger = logging.getLogger(__name__)


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


def make_config(default_map: Dict[str, Any], cls_name: str = "Config") -> Any:
    """Create a configuration dataclass.

    Creates a dataclass from a dictionary, and nested dictionary is transformed into
    another dataclass, which is subsequently nested in the main dataclass.

    Parameters
    ----------
    default_map: Dict[str, Any]
        configuration data from the CLI
    cls_name: str
        name of the dataclass

    Returns
    -------
    Config
        a configuration dataclass
    """
    data = default_map.copy()
    for key, value in data.items():
        if isinstance(value, dict):
            data[key] = make_config(value, key.title())(**value)

    return make_dataclass(cls_name, data)(**data)
