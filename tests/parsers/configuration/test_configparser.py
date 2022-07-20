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

from mdstab.parsers.configuration.configparser import Config
from mdstab.parsers.configuration.configparser import configure

from ...datafile import TRAJFILES


class TestConfig:
    """Test Config class."""

    @pytest.fixture
    def context(self) -> Config:
        """Create an mock click.Context object.

        Returns
        -------
        Config
            a mock click.Context object
        """
        config = Config()
        config.analysis = "coordinates"
        return config

    def test_config(self, context: Config) -> None:
        """Test the Config class.

        GIVEN a Config class
        WHEN the class is initialized
        THEN the object should be of type `dataclass`

        Parameters
        ----------
        context: Config
            mock object
        """
        assert is_dataclass(context)
        assert hasattr(context, "analysis")
        assert context.analysis == "coordinates"

    def test_update(self, context: Config) -> None:
        """Test functionality of update method.

        GIVEN a Config object
        WHEN a dict is included in the update method
        THEN the object should add the new attribute

        Parameters
        ----------
        context: Config
            mock object
        """
        kwargs = dict(debug=True, startres=1)
        context.update(**kwargs)

        assert hasattr(context, "analysis")
        assert hasattr(context, "debug")
        assert hasattr(context, "startres")
        assert context.debug
        assert context.startres == 1


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
