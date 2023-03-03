#!/usr/bin python3

""" Cancer Genome Atlas Parser/Mutation table """

"""
Program:    mutation_parser
File:       mutation_parser.py
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
This program parses gdc files (https://portal.gdc.cancer.gov) for kinesin database population.

Usage:
======
gdc_parser         SELF

Revision History:
=================
V1.2    22.01.19        PyMysql insert data, create dictionary      By: JJS
V1.3    23.01.19        Rewrite into functions                          JJS     
V1.4    23.01.19        Added functions for source and impact           JJS
                        tables
V1.5    25.01.19        Added file opening and gene specification to    JJS
                        parseMutation function. Removed insert          
                        functions to a new module. Renamed module.
"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import config_kinesin
import json

# ******************************************************************************
def parseMutations(gene, json_file):
    """ Function to parse json files from Genomic Data Commons.
    Input                   file                file of JSON mutations
    Output                  mut_list  [0]       list of strings for mutation table
                            source_list  [1]    list of strings for source table
                            impact_list  [2]    list of strings for impact table
    """

    mut_list    = []
    mut_entry   = []
    source_entry= []
    source_list = []
    impact_entry = []
    impact_list = []

    source_db   = 'GDC'
    my_gene     = 'KIF11'

    with open ('mutations.2018-10-03.json', 'r') as f:
        mutations = json.load(f)
        ## first index is the item number, 'x'
        x = 0
        for k in mutations:
            try:
                gene_name       = str(mutations[x]['consequence'][0]['transcript']['gene']['symbol'])
                protein         = str(mutations[x]['consequence'][0]['transcript']['aa_change'])
                consequence     = str(mutations[x]['consequence'][0]['transcript']['consequence_type'])
                genomic_id      = str(mutations[x]['genomic_dna_change'])
                ## eliminate 'g' in genomic_id
                genomic_id      = genomic_id.replace('g.', '')
                source_id       = str(mutations[x]['ssm_id'])
                mutation_type   = str(mutations[x]['mutation_subtype'])
                vep = str(mutations[x]['consequence'][0]['transcript']['annotation']['vep_impact'])
                if vep == '':
                    vep = ' '
                polyphen = str(mutations[x]['consequence'][0]['transcript']['annotation']['polyphen_impact'])
                if polyphen == '':
                    polyphen = ' '
                sift = str(mutations[x]['consequence'][0]['transcript']['annotation']['sift_impact'])
                if sift == '':
                    sift = ' '
                coding          = 'y'
                cds             = ' '
                organism        = 'Homo sapiens'
                domain          = ' '              # determine later
                x += 1
            except Error as e:
                print("Error", e)
            # put all strings into lists for each table
            if gene_name == my_gene:            ## check that gene is KIF11
                mut_entry       = [genomic_id, coding, cds, mutation_type, consequence, protein, gene_name, organism, domain]
                source_entry    = [source_id, source_db, protein]
                impact_entry    = [protein, vep, sift, polyphen]
                impact_list.append(impact_entry)
                mut_list.append(mut_entry)
                source_list.append(source_entry)
    f.close()

    return mut_list, source_list, impact_list

# ******************************************************************************

########## main ############

if __name__ == "__main__":

    gdc_att    = parseMutations('KIF11', 'mutations.2018-10-03.json')
    gdc_mut    = gdc_att[0]
    gdc_source = gdc_att[1]
    gdc_impact = gdc_att[2]

    print(gdc_att[0])

    # insertMutation(gdc_mut)
    # insertSource(gdc_source)
    # insertImpact(gdc_impact)

    i=0
    for k in gdc_mut:
        if k[5] != "None": # and k[4] == "missense_variant":
            i += 1
        #print(k)

    print('The total number of mutations is: ', i )






