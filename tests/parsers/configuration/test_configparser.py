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
"""Test ConfigParser class."""
from dataclasses import is_dataclass

import click
import pytest

from mdstab.parsers.configuration.configparser import configure
from mdstab.parsers.configuration.configparser import make_config

from ...datafile import TRAJFILES


class TestConfigParser:
    """Test ConfigParser class."""

    @pytest.fixture
    def ctx(self) -> click.Context:
        """Create a mock click.Context object.

        Returns
        -------
        Context
            a mock click.Context object
        """
        ctx = click.Context(click.Command("mock"))
        ctx.default_map = {}
        return ctx

    def test_configure(self, ctx: click.Context) -> None:
        """Test the configure function.

        GIVEN a click.Context and a filename
        WHEN the configure function is called
        THEN context.default_map is populated with the YAML configuration

        Parameters
        ----------
        ctx: click.Context
            mock object
        """
        configure(ctx, [], TRAJFILES)

        assert len(ctx.default_map) > 0  # type: ignore
        assert "trajfiles" in ctx.default_map  # type: ignore
        assert ctx.default_map["trajfiles"] is not None  # type: ignore

    def test_make_config(self, ctx: click.Context) -> None:
        """Test the make_config function.

        Parameters
        ----------
        ctx: click.Context
            mock object
        """
        configure(ctx, [], TRAJFILES)
        config = make_config(ctx.default_map)

        assert is_dataclass(config)
        assert hasattr(config, "trajfiles")
        assert config.trajfiles is not None
