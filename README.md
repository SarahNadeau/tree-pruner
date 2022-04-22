# Tree pruner

Select a small, representative set of sequences off of a larger phylogenetic tree.

This code selects sequences to represent the major clades of a phylogenetic tree.
The aim is that a query isolate could be placed on a phylogeny of these representative sequences for a quick first-look at the putative origins of the isolate.

## Requirements
* `ape` 
* `argparse` 
* `dplyr` 
* `phytools` 
* `tidytree`

## Usage
`prune_tree.R [-h] [-t] [-f TREEFILE] [-n NSEQS] [-o OUTDIR]`

### Options
```
  -h, --help            show this help message and exit
  -t, --test            do a test run
  -f TREEFILE, --treefile TREEFILE
                        tree to choose representative sequences from
  -n NSEQS, --nseqs NSEQS
                        number of sequences to select [default 100]
  -o OUTDIR, --outdir OUTDIR
                        directory for results [default "results"]
```
