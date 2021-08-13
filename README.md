# hgfind
[![pytest](https://github.com/mnahinkhan/hgfind/actions/workflows/python-package.yml/badge.svg)](https://github.com/mnahinkhan/hgfind/actions/workflows/python-package.yml)
[![black](https://github.com/mnahinkhan/hgfind/actions/workflows/black-check.yml/badge.svg)](https://github.com/mnahinkhan/hgfind/actions/workflows/black-check.yml)
[![isort](https://github.com/mnahinkhan/hgfind/actions/workflows/isort-check.yml/badge.svg)](https://github.com/mnahinkhan/hgfind/actions/workflows/isort-check.yml)
[![PyPI version](https://badge.fury.io/py/hgfind.svg)](https://badge.fury.io/py/hgfind)

Get the human genome (hg38) coordinates of a gene

## Installation
```
pip install hgfind
```

*Requirements*:
 - Python >=3.6

## Usage

`hgfind` can be used either on the command line or as a function within
Python

### As a commandline tool
```
hgfind <gene>
```

where `gene` is a gene name such as "PTEN". The name can be a synonym as well.

A successful example:
```
$ hgfind auf1
HNRNPD => 4:82352498-82374503 (-)

$ echo $?
0
```

The result shows that `HNRNPD` (a synonym for `AUF1`) lies on chromosome 4,
in the specified base interval. The `(-)` indicates that its transcripts
all lie on the reverse strand.


Using an unrecognized name results in an error:
```
$ hgfind fjlsfl
fjlsfl not recognized as a gene

$ echo $?
1
```

### As a function in Python
As an example on the Python REPL:
```
>>> from hgfind import hgfind
>>> hgfind("Neat2")
{'chr_n': 11, 'start_coord': 65497688, 'end_coord': 65506516, 'strand': '+', 'official_name': 'MALAT1'}
```

Using an unrecognized name results in an error:
```
>>> hgfind("gewgwre")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "hgfind.py", line 140, in hgfind
    raise WrongGeneName(
hgfind.WrongGeneName: {'message': 'The input gene could not be recognized', 'gene': 'GEWGWRE'}
```


## Contributing
Any suggestions / PR requests are welcome!

## Development
Enable recommended Git Hooks as follows:
```
git config --local core.hooksPath .githooks/
```
The above will run the following to ensure code consistency every time you
commit:
 - [black](https://github.com/psf/black)
 - [isort](https://github.com/PyCQA/isort)

Also use [fit-commit](https://github.com/m1foley/fit-commit) to ensure
consistent commit message style.
