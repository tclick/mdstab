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
"""Various data files for testing."""
from pathlib import Path

from pkg_resources import resource_filename


__all__ = ["PROJ", "PROJNP", "TOPWW", "TRJWW"]

PROJ = resource_filename(__name__, Path().joinpath("data", "projection.csv").as_posix())
PROJNP = resource_filename(
    __name__, Path().joinpath("data", "projection.npy").as_posix()
)
TOP = resource_filename(__name__, Path().joinpath("data", "1qmt.prmtop").as_posix())
TOPWW = resource_filename(__name__, Path().joinpath("data", "protein.parm7").as_posix())
TRJ = resource_filename(__name__, Path().joinpath("data", "prod.nc").as_posix())
TRJWW = resource_filename(__name__, Path().joinpath("data", "protein.nc").as_posix())

# Cluster data
CENTROID = resource_filename(
    __name__, Path().joinpath("data", "pca-centroid.csv").as_posix()
)
CENTNPY = resource_filename(
    __name__, Path().joinpath("data", "pca-centroid.npy").as_posix()
)
CLUSTER = resource_filename(
    __name__, Path().joinpath("data", "pca-cluster.csv").as_posix()
)
CLUSTNPY = resource_filename(
    __name__, Path().joinpath("data", "pca-cluster.npy").as_posix()
)
LABELS = resource_filename(
    __name__, Path().joinpath("data", "pca-labels.npy").as_posix()
)
FRAMES = resource_filename(__name__, Path().joinpath("data", "frames.csv").as_posix())

TRAJFORM = resource_filename(
    __name__, Path().joinpath("data", "config-trajform.yaml").as_posix()
)

TRAJFILES = resource_filename(
    __name__, Path().joinpath("data", "config-trajfiles.yaml").as_posix()
)

TRAJDIH = resource_filename(
    __name__, Path().joinpath("data", "config-dihedrals.yaml").as_posix()
)

BAD_CONFIG = resource_filename(
    __name__, Path().joinpath("data", "config-bad.yaml").as_posix()
)
