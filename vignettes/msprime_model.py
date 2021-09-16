#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 21:23:26 2020
Modified on Mon Sep 7 16:51:08 2020

@author: alanaw

This script generates arrays that have dependent columns 
and exchangeable rows
"""

import msprime
import numpy as np

# Generate a genotype matrix under panmixia
# Simulates about 40,000 loci with recombination across the 
# loci
# Number of segregating sites can be estimated with Watterson's 
# formula. 
# This generates about 200,000 SNPs along a chromosome of length 5e7.
def getHaplotypes(N, P):
    sample_size = int(N)
    chr_length = int(P)  # set to 5e7
    ts =  msprime.simulate(Ne=1e4, 
        mutation_rate=2e-8,
        length=chr_length,
        recombination_rate=2e-8,
        sample_size=sample_size, 
        num_replicates=1)
    gts = []
    for tree_seq in ts:
        gts.append(tree_seq.genotype_matrix())
        gts = np.transpose(np.vstack(gts).astype(float))
    return(gts)


#print(gts)