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
import logging.config
import os
from typing import TYPE_CHECKING
from typing import Union

import yaml
from click_extra.logging import logger as clog
from license.licenses import GPLv3LaterLicense


if TYPE_CHECKING:
    PathLike = Union[str, os.PathLike[str]]
else:
    PathLike = Union[str, os.PathLike]

__version__: str = "0.0.0"
_name: str = "Timothy H. Click, Ph.D."
_email: str = "Timothy.Click@briarcliff.edu"
__copyright__: str = GPLv3LaterLicense.header(name=_name, email=_email)


def create_logger(
    name: str, level: int = logging.INFO, logfile: PathLike = "output.log"
) -> clog:
    """Create a logger.

    Parameters
    ----------
    name : str
        name of module
    level : int
        logging level
    logfile : PathLike
        filename

    Returns
    -------
    logger

    Raises
    ------
    ValueError
        If a filename is not defined.
    """
    if not logfile:
        raise ValueError("Filename not defined.")

    logger_dict = f"""
    version: 1
    disable_existing_loggers: False
    formatters:
        default:
            format: '%(asctime)s %(levelname)-8s %(name)-15s %(message)s'
            datefmt: '%Y-%m-%d %H:%M:%S'
    handlers:
        console:
            class: logging.StreamHandler
            level: {level}
            stream: ext://sys.stdout
            formatter: default
        file:
            class: logging.FileHandler
            filename: {logfile}
            level: {level}
            mode: w
            formatter: default
    root:
        level: {level}
        handlers:
            - console
            - file
    """
    logdict = yaml.safe_load(logger_dict)
    logging.config.dictConfig(logdict)
    return logging.getLogger(name)
