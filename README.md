# Molecular Dynamics Trajectory Analysis

[![PyPI](https://img.shields.io/pypi/v/mdstab.svg)][pypi_]
[![Status](https://img.shields.io/pypi/status/mdstab.svg)][status]
[![Python Version](https://img.shields.io/pypi/pyversions/mdstab)][python version]
[![License](https://img.shields.io/pypi/l/mdstab)][license]

[![Read the documentation at https://mdstab.readthedocs.io/](https://img.shields.io/readthedocs/mdstab/latest.svg?label=Read%20the%20Docs)][read the docs]
[![Tests](https://github.com/tclick/mdstab/workflows/Tests/badge.svg)][tests]
[![Codecov](https://codecov.io/gh/tclick/mdstab/branch/main/graph/badge.svg)][codecov]

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)][pre-commit]
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)][black]

[pypi_]: https://pypi.org/project/mdstab/
[status]: https://pypi.org/project/mdstab/
[python version]: https://pypi.org/project/mdstab
[read the docs]: https://mdstab.readthedocs.io/
[tests]: https://github.com/tclick/mdstab/actions?workflow=Tests
[codecov]: https://app.codecov.io/gh/tclick/mdstab
[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black

## Features

- Calculate r.m.s.d. of a protein
- Calculate r.m.s.f. of a protein
- Perform spatial or spatio-temporal decorrelation
  - Aligned Cartesian coordinates
  - `φ/ψ` dihedral angles
- Perform fluctuation matching
  - Models for various biomolecules
  - Cluster analysis
  - Spatial or spatio-temporal decorrelation

## Requirements

- Python 3.8+ (https://www.python.org)
- MDAnalysis 2.1+ (https://www.mdanalysis.org)
- modin (https://github.com/modin-project/modin)
- Academic charmm (https://www.academiccharmm.org)

## Installation

You can install _Molecular Dynamics Trajectory Analysis_ via [pip] from [PyPI]:

```console
$ pip install mdta
```

## Usage

Please see the [Command-line Reference] for details.

## Contributing

Contributions are very welcome.
To learn more, see the [Contributor Guide].

## License

Distributed under the terms of the [GPL 3.0 license][license],
_Molecular Dynamics Trajectory Analysis_ is free and open source software.

## Issues

If you encounter any problems,
please [file an issue] along with a detailed description.

## Credits

This project was generated from [@cjolowicz]'s [Hypermodern Python Cookiecutter] template.

[@cjolowicz]: https://github.com/cjolowicz
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/tclick/mdstab/issues
[pip]: https://pip.pypa.io/

<!-- github-only -->

[license]: https://github.com/tclick/mdstab/blob/main/LICENSE
[contributor guide]: https://github.com/tclick/mdstab/blob/main/CONTRIBUTING.md
[command-line reference]: https://mdstab.readthedocs.io/en/latest/usage.html
