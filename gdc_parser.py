#!/usr/bin python3

""" Cancer Genome Atlas Parser """

"""
Program:    gdc_parser
File:       gdc_parser.py
Version:    1.0
Date:       03.10.18
Function:   Parse gdc files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses gdc files for kinesin database population

Usage:
======
gdc_parser         SELF

Revision History:
=================


"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import json
from pprint import pprint

# ******************************************************************************

# parse genomic data commons files

with open ('mutations.2018-10-03.json', 'r') as f:
    mutations = json.load(f)

#pprint(mutations)

## first index is the item number
## check that gene is KIF11
mut_dict = {}
my_gene = 'KIF11'
x = 0
for y in mutations:
    try:
        gene        = mutations[x]['consequence'][0]['transcript']['gene']['symbol']
        aa_change   = mutations[x]['consequence'][0]['transcript']['aa_change']
        cons_type   = mutations[x]['consequence'][0]['transcript']['consequence_type']
        vep         = mutations[x]['consequence'][0]['transcript']['annotation']['vep_impact']
        genomic     = mutations[x]['genomic_dna_change']
        mut_type    = mutations[x]['mutation_subtype']
    except IndexError:
        break
    x += 1
    if gene == my_gene:
        ## iterate thru json file and insert data for each mutation entry
        ## pymysql script 
        #print(genomic, aa_change, cons_type)
        ## could make a dictionary of mutations and attributes. too many?
        mut_dict[genomic] = (aa_change, cons_type, mut_type)       
    else:
        print('Gene not found!')
        sys.exit()
print('The total number of mutations is: ', x)
for x in mut_dict:
    print(x, mut_dict[x])
