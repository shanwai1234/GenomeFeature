# GenomeFeature [![](https://img.shields.io/badge/Release-v1.0.1-blue.svg)](https://github.com/shanwai1234/GPWAS/commits/master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
-------------

GenomeFeature is a python-based package which is used for extracting genome and protein features of each single annotated gene, and is used for Machine Learning/Deep Learning purpose.

## Dependencies

- [BioPython](https://biopython.org/wiki/Download)
- [pyfasta](https://pypi.org/project/pyfasta/)

## How to use

1. Download this python-based package to local
2. Entering the downloaded package by
```
$ cd GenomeFeature
```
3. Install this package to local environment by typing
```
$ sudo pip install -e .
```
4. Preparing 4 required files for later feature generation, which are cds sequences of primary transcript, assembled genome, annotation file and amino acid sequences of primary transcript. Then moving all of four files into one folder.
5. To use this package, either seeing GenomeFeature/GenomeFeature/genomefeature.py or using main functions in this package.
