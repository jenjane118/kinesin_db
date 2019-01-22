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
22.01.19        PyMysql insert data, create dictionary      JJS
     

"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import config_kinesin
import json
from pprint import pprint

# ******************************************************************************

# parse genomic data commons files

#Connect to MySQL database
cnx = pymysql.connect(host  =config_kinesin.database_config['dbhost'],
                      port  =config_kinesin.database_config['port'],
                      user  =config_kinesin.database_config['dbuser'],
                      passwd=config_kinesin.database_config['dbpass'],
                      db    =config_kinesin.database_config['dbname'])
cursor = cnx.cursor(pymysql.cursors.DictCursor)



with open ('mutations.2018-10-03.json', 'r') as f:
    mutations = json.load(f)


## first index is the item number
## check that gene is KIF11
mut_list = []
mutation_dict ={}
my_gene = 'KIF11'
x = 0
for y in mutations:
    try:
        gene_name       = str(mutations[x]['consequence'][0]['transcript']['gene']['symbol'])
        protein         = str(mutations[x]['consequence'][0]['transcript']['aa_change'])
        consequence     = str(mutations[x]['consequence'][0]['transcript']['consequence_type'])
 #       vep            = mutations[x]['consequence'][0]['transcript']['annotation']['vep_impact']
        genomic_id      = str(mutations[x]['genomic_dna_change'])
        ## eliminate 'g' in genomic_id
        genomic_id      = genomic_id.replace('g.', '')
        mutation_type   = str(mutations[x]['mutation_subtype'])
    except IndexError:
        break
    if gene_name == my_gene:
        x = x+1
        mut_list = [genomic_id, 'y', ' ', mutation_type, consequence, protein, gene_name, 'homo sapiens']
        if mut_list[0] != "None" and mut_list[4]=="missense_variant":
            # create dictionary of missense mutations
            mutation_dict[mut_list[5]] = mut_list[0:5]
            # iterate thru json file and insert data for each missense mutation entry
            sql_mutation = "INSERT INTO mutation (genomic, coding, cds, mutation_type, consequence, protein, " \
                "gene_name, organism) VALUES(%s,%s,%s,%s,%s,%s,%s,%s)"
            rows = cursor.execute(sql_mutation, mut_list)

    else:
        print('Gene not found')
        sys.exit()

cnx.commit()
cnx.close()

i = 0
for h in sorted(mutation_dict):
    i += 1
    print(h, mutation_dict[h])
print('The total number of missense mutations is: ', i )


# next eliminate 'chr' before genomic_id, find info for cds column?
# remember to delete all in mutation table on mysql database before running again