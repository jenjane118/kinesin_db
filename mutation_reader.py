#!/usr/bin python3

""" GDC API module """

"""
Program:    gdc_api_module
File:       gdc_api_module.py
Version:    1.0
Date:       04.03.19
Function:   Utilise TGDC API to get more info for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses gdc files tsv files (https://portal.gdc.cancer.gov) for use in API.

Usage:
======
id_reader                SELF

Revision History:
=================
V 1.0                    Initial                    By: JJS

"""

# ******************************************************************************
# Import libraries
import sys
import csv
import json

# ******************************************************************************

def IdReader(csv_file):
    """This function will read the mutation ids from downloaded mutations from TGDC assoc w/ KIF11 gene
    Input               tsv_file                    Downloaded from website
    Output              id_list                     List of mutation ids for submission to API (get or post?)
    """

    mutation_dict = {}

    # open csv file
    with open(csv_file, 'r') as file:
        csv_file = csv.reader(file, delimiter=' ')
        for row in csv_file:
            mutation_dict['ids'] = row


    with open('data.json', 'w') as outfile:
                json.dump(mutation_dict, outfile, sort_keys = True, indent = 4)


    return outfile

# *************************************************************************************

# def IdsJson(tsv_file):
#     with open(tsv_file) as file:
#




# *************************************************************************************
########## main ############

if __name__ == "__main__":

    tgdc = IdReader('tgdc_mutationset.tsv')

    print(tgdc)