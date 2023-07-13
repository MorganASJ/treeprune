#!/bin/bash
iqtree2 -s *.pruned.fas -m GTR+F+G -z all_trees.pruned.treels -te *02_concat_tree.phy.contree.pruned.tre -zb 10000 -au

exit 0
