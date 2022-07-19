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
"""Molecular Dynamics Setup and Trajectory Analysis for Biomolecules."""
import logging
import os
from typing import TYPE_CHECKING
from typing import Any
from typing import Dict
from typing import Union


logger: logging.Logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

if TYPE_CHECKING:
    PathLike = Union[str, os.PathLike[str]]
else:
    PathLike = Union[str, os.PathLike]

__version__: str = "0.0.0"

_MASK: Dict[str, str] = dict(
    ca="protein and name CA",
    cab="protein and name C[AB]",
    backbone="backbone",
    sidechain="protein and not backbone and not element H*)",
    noh="protein and not element H*)",
    all="all",
)


def create_logging_dict(logfile: PathLike, level: int) -> Dict[str, Any]:
    """Configure the logger.

    Parameters
    ----------
    logfile : PathLike
        Filename for log output.
    level : int
        logging level
    Returns
    -------
    Dict
        Configuration data for logging.

    Raises
    ------
    ValueError
        If a filename is not defined.
    """
    if not logfile:
        raise ValueError("Filename not defined.")

    logger_dict = dict(
        version=1,
        disable_existing_loggers=False,  # this fixes the problem
        formatters=dict(
            standard={
                "class": "logging.Formatter",
                "format": "%(name)-12s %(levelname)-8s %(message)s",
            },
            detailed={
                "class": "logging.Formatter",
                "format": ("%(asctime)s %(name)-15s %(levelname)-8s " "%(message)s"),
                "datefmt": "%m-%d-%y %H:%M",
            },
        ),
        handlers=dict(
            console={
                "class": "logging.StreamHandler",
                "level": level,
                "formatter": "standard",
            },
            file={
                "class": "logging.FileHandler",
                "filename": logfile,
                "level": level,
                "mode": "w",
                "formatter": "detailed",
            },
        ),
        root=dict(level="INFO", handlers=["console", "file"]),
    )
    return logger_dict
