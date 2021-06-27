# hgfind

Get the human genome (hg38) coordinates of a gene

## Installation
```
pip install hgfind
```

## Usage

`hgfind` can be used either on the commandline or as a module within
Python

### As a commandline tool
```
hgfind <gene>
```

where `gene` is a gene name such as "PTEN". The name can be a synonym as well.

A successful example:
```
$ hgfind auf1
HNRNPD => 4:82352498-82374503

$ echo $?
0
```

Using an unrecognized name results in an error:
```
$ hgfind fjlsfl
fjlsfl not recognized as a gene

$ echo $?
1
```



### As a Python module
As an example on the Python REPL:
```
>>> from hgfind import hgfind
>>> hgfind("Neat2")
{'chr_n': 11, 'start_coord': 65497688, 'end_coord': 65506516, 'official_name': 'MALAT1'}
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