#!/usr/bin python3

""" Cancer Genome Atlas Parser/Mutation table """

"""
Program:    cosmic_parser
File:       cosmic_parser.py
Version:    1.0
Date:       03.10.18
Function:   Parse cosmic files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses COSMIC files (https://cancer.sanger.ac.uk/cosmic/).

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
#import pymysql
#import config_kinesin
import csv
import re

# ******************************************************************************

def cosmicParser(file):
    """Function to parse mutation files from COSMIC database
    Input               file                    csv file download
    Output              id_list                 list of mutation_id's
                        mutation_dict           dictionary of attributes of mutation and impact
                        source_list             cosmic mutation entries for source table

    """
    mut_list        = []
    mut_entry       = []
    mutation_dict   = {}

    my_gene = 'KIF11'

    id_entry    = []
    id_list     = []

    # parse out info into attribute names
    x = 0
    for row in file:
        try:
            # mutation table  (only include info not available in gdc files)
            gene_name       = row[0]
            #genomic_id      = row[1]
            #coding          = 'y'
            protein         = row[18]
            # eliminate 'p.'
            protein         = protein.replace('p.', '')
            cds             = row[17]
            #mutation_type   = row[19]           # need all before hyphen
            #consequence     = row[19]           # after the hyphen
            #organism        = 'Homo sapiens'
            #domain          = ' '               # calculate later
            # source_info table
            source_db       = 'COSMIC'
            source_id       = row[16]
            # impact table
            fathhm_score    = row[28]
            fathhm_pred     = row[27]
            # tissue table
            tissue_type     = row[7]
            cancer_type     = row[12]
        except Error as e:
            print("Error", e)

        # make dictionary to fill in missing attribute fields in mutation table before insertion
        if protein != '?' and gene_name == my_gene:  ## check that gene is KIF11 don't include ambiguous amino acid changes
            mutation_dict[protein] = (cds, source_db, source_id, fathhm_score, fathhm_pred, tissue_type,
                                            cancer_type)

    return mutation_dict

# ******************************************************************************
##########  MAIN ####################

if __name__ == "__main__":

    with open('V87_38_MUTANT.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        # to determine index number of columns in csv file
        # line_count = 0
        # for row in csv_reader:
        #     i = 0
        #     if line_count == 0:
        #         for r in row:
        #             print(i, r)
        #             i +=1
        #         line_count += 1
        #     else:
        #         line_count += 1
        cosmic_dict = cosmicParser(csv_reader)
        for x in cosmic_dict:
            print(x, cosmic_dict[x])


# ******************************************************************************
