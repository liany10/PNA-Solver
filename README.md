![logo](logo.png)

## *PNA Solver* 

#### - Python code for solving PNA issues


## Contact

#### LIAN Yun
#### Chinese University of Hong Kong, (Shen Zhen)
#### Email: 219059026@link.cuhk.edu.cn

## Requirement: 
#### Environment: Python3


#### Dependencies:

| Dependencies | Version >= |
| ------------- |:-------------:|
| re | 3.1.2 |
| subprocess | 3.1.2 |


## Introduction

PNA solver has two modules, sequence.py and PNA concentration calculator.py. The sequence.py provides some basic functions on sequence operations for DNA, RNA and PNA. The PNA concentration calculator.py provides some more complicated functions dealing with data bank including searching for possible RNA targets according to PNA's structure. It involves RNAfold from the Vienna RNA package.

By importing PNA concentration calculator.py, one can accomplish tasks including:
1) Calculating extinction coefficient and molecular weight by XNADealer.ReadPNAandCalculate().
2) Searching for exact complementary sequence in RNAs for PNAs by XNADealer.FindExactPairsInSpecies().
3) Searching for RNAs that have certain secondary structures for PNA targeting by XNADealer.FindRNAsFordbPNA().


## Reference

Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
ViennaRNA Package 2.0
Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
