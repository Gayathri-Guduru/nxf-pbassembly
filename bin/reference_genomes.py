#!/usr/bin/env python

"""
This program takes the sample species list as input 
and creates table with reference genome .fasta and .gff files
required for assembly QC
"""

import os.path



os.path.join(
    's3://zymo-sequencing/r64524e_$(date)_*/1_A01/', 'in', 'here')
