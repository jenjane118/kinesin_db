#!/usr/bin python3

""" Clinvar Parser Module """

"""
Program:    clinvar_impact
File:       clinvar_impact.py
Version:    1.0
Date:       10.04.19
Function:   Parse Clinvar files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses Clinvar file (https://www.ncbi.nlm.nih.gov/clinvar) for clinical significance. This file includes
all recorded mutations for the gene 'KIF11' that involve only a single gene. Most of these are germline variants.
It includes functions to parse the output file from clinvar and update mysql database tables with clinical significance.

Usage:
======
vep_parse                SELF

Revision History:
=================
V1.0    10.04.19        Initial                                     By: JJS

"""

#
#*****************************************************************************
# Import libraries

import csv
import re
import config_home
import pymysql

#*****************************************************************************

def parseClinvar(gene, csv_file):
    """This function parses Clinvar results file.
    Input               gene                KIF11
                        csv_file            file of Clinvar results returned from web query
    Output              cv_list             list of mutation identifier and clinical significance
    """

    cv_entry = []
    cv_list  = []
    # open csv file
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        first_row = True
        for row in csv_reader:
            if first_row == True:           #skip first header row with titles
                first_row = False
            else:
                try:
                    gene_name   = row[1]
                    genomic     = row[8]    #genomic location on chromosome 10
                    # remove anything after first position number
                    p       = re.compile(r'\s\-.*$')
                    genomic = p.sub('', genomic)
                    # reformat genomic entry to include chromosome number
                    genomic     = '10:' + genomic
                    clin_sig    = row[3]    #clinical significance
                    # remove review date in parentheses
                    q           = re.compile(r'\(.*\)$')
                    clin_sig    = q.sub('', clin_sig)

                except Exception as e:
                    print("Error", e)

                if gene_name == gene:       #check that gene is correct one
                    cv_entry = [genomic, clin_sig]
                    cv_list.append(cv_entry)

    return cv_list

# ******************************************************************************

def clinvarUpdate(clinvar_list):
    """ This function uses list of clinvar attributes to update impact table in kinesin database.
    Input               clinvar_list                    list of attributes (genomic location, clinical significance)
    Output              count                           number of successfully updated rows in impact table
    """

    # Connect to MySQL Database
    cnx = pymysql.connect(host=config_home.database_config['dbhost'],
                          port=config_home.database_config['port'],
                          user=config_home.database_config['dbuser'],
                          passwd=config_home.database_config['dbpass'],
                          db=config_home.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    impact_list = []
    # use genomic attribute to find protein(mutation_id) for impact table
    query = "SELECT protein FROM mutation WHERE genomic = %s ;"
    for x in clinvar_list:
        mutation = str.format(x[0])         # retrieve mutation (genomic location) from parsed clinvar list
        with cnx.cursor() as cursor:
            cursor.execute(query, (mutation))
            temp = cursor.fetchone()

        x.pop(0)        # remove genomic location from first position
        x.append(temp)  # add mutation_id (aa change) to last position
        impact_list.append(x)

    sql_clinvar = "UPDATE impact SET clinvar_prediction = %s WHERE mutation_id = %s;"

    count = 0
    for item in impact_list:
        with cnx.cursor() as cursor:
            rows = cursor.execute(sql_clinvar, item)
            count += 1

    cnx.commit()
    cnx.close()

    return count

# ******************************************************************************

########## main ############

if __name__ == "__main__":

    mygene = 'KIF11'
    myfile = 'clinvar_result.txt'

    results = parseClinvar(mygene, myfile)
    #print(results)

    successful = clinvarUpdate(results)
    print(successful)