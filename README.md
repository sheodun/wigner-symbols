# Wigner-Symbols
A simple zero-dependenct module for calculating Wigner 3j and 6j symbols.

## Installation

### Install package from pip

This package has not yet been published on PyPI.

### Local installation

Navigate to the directory and run

```bash
pip install .
```

## Usage

It is very simple to calculate 3j and 6j symbols:

```python
from wigner_symbols import calculate_3j, calculate_6j

calculate_3j(2.0, 2.0, 2.0, 0.0, 0.0, 0.0)  # returns -0.2390457
calculate_6j(1.0, 2.0, 3.0, 2.0, 1.0, 2.0)  # returns 0.04364358
```

## Development setup

Install all dependencies:

```bash
poetry install
```

To run tests:

```bash
poetry run pytest
```

## Meta

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/sheodun](https://github.com/sheodun)

## Contributing

1. Fork this repository
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
