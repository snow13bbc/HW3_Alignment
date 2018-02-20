# HW3_Alignment
Homework 3 Assignment: 

[![Build
Status](https://travis-ci.org/snow13bbc/HW3_Alignment.svg?branch=master)](https://travis-ci.org/snow13bbc/HW3_Alignment)

Travis - https://travis-ci.org/snow13bbc/HW3_Alignment

Skeleton for alignment optimization project.

## assignment

1. Implement the Smith-Waterman algorithm. 
2. Test different gap penalties and scoring matrices.
3. Optimize algorithm. 



## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `HW3_Alignment/__main__.py`) can be run as
follows

```
python -m Alignment 
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.
