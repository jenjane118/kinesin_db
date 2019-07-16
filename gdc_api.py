#!/usr/bin python3

""" GDC API module """

"""
Program:    gdc_api
File:       gdc_api.py
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
This program parses gdc files (https://portal.gdc.cancer.gov) for 'ssmid' and uses 
to make API requests for tissue and cancer type information.

Usage:
======
id_reader                SELF

Revision History:
=================
V 1.0        06.03.19            Initial                     By: JJS
V 1.1        07.03.19            Add request script              JJS
"""

# ******************************************************************************
# Import libraries
import sys
import csv
import json
import requests

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

def GdcRequest(file, results_file):
    """This function will submit a request for each of the mutation id's and retrieve tissue info from webservice
    Input               file                    List of mutational ids related to gene (json)
                        results_file            Name of results file
    Output              outfile                 Dictionary (json) with response data
    """

    tissue_list = []

    cases_endpt = 'https://api.gdc.cancer.gov/ssms/'
    moreparams = '?expand=occurrence.case&pretty=true'
    # testids = ["8604d5d0-fd49-5452-935a-f69a6c24de46"]

    with open(file) as jsons:
        configs = json.load(jsons)

    for x in configs:
        for i in configs[x]:
            if i != 'id':
                for_request = ''.join([cases_endpt, i, moreparams])
                r = requests.get(for_request)
                rj = r.json()
                tissue_list.append(rj)
    with open(results_file, 'w') as outfile:
        json.dump(tissue_list, outfile, sort_keys=True, indent=4)

    return outfile

# *************************************************************************************
########## main ############

if __name__ == "__main__":

    tgdc = IdReader('tgdc_mutationset.tsv')

    tissue_data = GdcRequest('data.json', 'results.json')
