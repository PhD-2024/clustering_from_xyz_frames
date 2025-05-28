## How to use
# What does it do?
This script takes an xyz file and uses a simple cutoff criterion to get "molecules" or "fragments" of more than 1 atoms.
Single atoms are most likely metal centers and get ignored.
For the future also pbc structures may be supported. If desired tell me and I will add a pbc treatment for the distance
determination.

# Requirements

python3 - nothing else 

# How to run
example:

```
python3 read_and_cluster.py --filename testfile_h.xyz --cutoff 1.65 --indexing_1 False --cutoff 1.65
```

# Argparser options

filename: your input file
cutoff: the cutoff you want to use
indexing_1: whether you want 1 or zero indexing (False for 0)
outname: string - outputs with more than 1 atom will be returned in outname_(NR).txt
