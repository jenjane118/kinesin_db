#!/usr/bin python3

""" Cancer Genome Atlas Parser/Mutation table """

"""
Program:    cosmic_parser
File:       cosmic_parser.py
Version:    1.0
Date:       23.01.19
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
V 1.0               Initial version                                 By: JJS
V 1.2               Includes file opening and gene specification        JJS
                    as function arguments.            
                      
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

def cosmicParser(my_gene, csv_file):
    """Function to parse mutation files from csv files from COSMIC database
    Input               file                    csv file download
    Output              id_list                 list of mutation_id's
                        mutation_dict           dictionary of attributes of mutation and impact

    """

    mut_entry       = []
    mutation_dict   = {}

    # open csv file
    with open(csv_file) as file:
        csv_reader = csv.reader(file, delimiter=',')
        # parse out info into attribute names
        for row in csv_reader:
            try:
                # mutation table  (only include info not available in gdc files)
                gene_name       = row[0]
                genomic_id      = row[1]
                coding          = 'y'
                protein         = row[18]
                # eliminate 'p.'
                protein         = protein.replace('p.', '')
                cds             = row[17]
                description     = row[19]             # includes type and consequence
                organism        = 'Homo sapiens'
                domain          = ' '                 # calculate later

                # source_info table
                source_db       = 'COSMIC'
                source_id       = row[16]
                # impact table
                fathhm_score    = row[28]
                fathhm_pred     = row[27]

                # tissue table
                tissue_type     = row[7]
                cancer_type     = row[12]
                if cancer_type  == 'NS':
                    cancer_type = row[11]
            except Error as e:
                print("Error", e)

            p   = re.compile(r'(^\w+)\s-\s(.+$)')
            it  = p.finditer(description)
            for match in it:
                mutation_type   = match.group(1)  # need all before hyphen
                consequence     = match.group(2)  # after the hyphen

            # make dictionary to fill in missing attribute fields in mutation table before insertion
            if protein != '?' and gene_name == my_gene:  ## check that gene is KIF11, don't include ? for aa_change
                mutation_dict[protein] = (cds, mutation_type, consequence)

    return mutation_dict

# ******************************************************************************
##########  MAIN ####################

if __name__ == "__main__":

    cosmic_dict = cosmicParser('KIF11', 'V87_38_MUTANT.csv')
    for x in cosmic_dict:
        print(x, cosmic_dict[x])


# ******************************************************************************
