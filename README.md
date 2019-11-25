# GenomeFeature [![](https://img.shields.io/badge/Release-v1.0.2-blue.svg)](https://github.com/shanwai1234/GPWAS/commits/master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
-------------

GenomeFeature is a python-based package which is used for extracting genome and protein features of each single annotated gene, and is used for Machine Learning/Deep Learning purpose.

## Dependencies

- [BioPython](https://biopython.org/wiki/Download)
- [pyfasta](https://pypi.org/project/pyfasta/)

## How to use

1. Download this python-based package to local; Enter the `GenomeFeature` folder
2. Build and install
```
$ python setup.py install
```
3. Install this package to local environment by typing
```
$ sudo pip install -e .
```
4. Enter any directory containing your gff and fasta file and run command line like below
```
$ python -m genomefeature -a Zmays_493_RefGen_V4.gene_exons_xm_primarytranscript_only.gff3 -g Zmays_493_APGv4.softmasked-fixname.fa -d 200 -u 200 -e 1000 -o maize-test.txt

-a gff file
-g fasta file
-d 3'-UTR length
-u 5'-UTR length
-e extension length from UTR regions in both upstream and downstream
-o output file of features 
```
Note: Applying this code to any genome needs a well-organized gff and fasta file. 
1> header in fasta file should be consistent with chromosome number in gff file; 
2> example in gff file should be in this following format (if your gene name contain '.', removing '.' sign in the gene name in the last column of gff file, this will cause error for parsing gene name):
```
1       phytozomev12    gene    44289   49837   .       +       .       ID=Zm00001d027230.RefGen_V4;Name=Zm00001d027230
1       phytozomev12    mRNA    44289   49837   .       +       .       ID=Zm00001d027230_T001.RefGen_V4;Name=Zm00001d027230_T001;pacid=40214039;longest=1;Parent=Zm00001d027230.RefGen_V4
1       phytozomev12    exon    44289   44947   .       +       .       ID=Zm00001d027230_T001.RefGen_V4.exon.1;Parent=Zm00001d027230_T001.RefGen_V4;pacid=40214039
1       phytozomev12    five_prime_UTR  44289   44350   .       +       .       ID=Zm00001d027230_T001.RefGen_V4.five_prime_UTR.1;Parent=Zm00001d027230_T001.RefGen_V4;pacid=40214039
1       phytozomev12    CDS     44351   44947   .       +       0       ID=Zm00001d027230_T001.RefGen_V4.CDS.1;Parent=Zm00001d027230_T001.RefGen_V4;pacid=40214039
```

