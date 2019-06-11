#!/usr/bin python3

""" FATHMM Parser Module """

"""
Program:    fathmm_impact
File:       fathmm_impact.py
Version:    1.0
Date:       14.05.19
Function:   Parse FATHMM files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses FATHMM results file downloaded from FATHMM website (http://fathmm.biocompute.org.uk) for predictions
of the oncogenic status of missense amino-acid substitutions in the KIF11 protein. This file includes
all missense mutations for the gene 'KIF11' and assesses 'driver' or 'passenger' status. The prediction threshold for 
the algorithm was set for '1.0', which is slightly more sensitive and less-specific than the default values.

Usage:
======
vep_parse                SELF

Revision History:
=================
V1.0    14.05.19        Initial                                     By: JJS

"""

#
#*****************************************************************************
# Import libraries

import csv
import re
import config_home
import pymysql
import config_kinesin

#*****************************************************************************

def fathmmResultsParser(csv_file):
    """This function parses FATHMM results file.
        Input               csv_file            file of FATHMM results returned from web query
        Output              fathmm_dict         dictionary mutation identifier:prediction, score
        """
    fathmm_list  = []
    fathmm_entry = []
    fathmm_dict  = {}
    with open (csv_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        row_one = True
        for row in csv_reader:
            if row_one == False:
                try:
                    mutation     = row[3]
                    prediction   = row[4]
                    score        = row[5]
                    fathmm_dict[mutation] = prediction, score
                except NameError as e:
                    print("Error", e)
            else:
                row_one = False
    file.close()

    return fathmm_dict

#*****************************************************************************

def fathmmInsert(results_dict, database):
    """ This function uses list of clinvar attributes to update impact table in kinesin database.
        Input               results_dict                    dict of attributes (mutation_id:prediction, score)
                            database                        which database used (home or kenobi)
        Output              count                           number of successfully updated rows in impact table
        """
    # Connect to MySQL Database (kinesin on kenobi)
    if database == 'kenobi':
        cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                              user=config_kinesin.database_config['dbuser'],
                              passwd=config_kinesin.database_config['dbpass'],
                              db=config_kinesin.database_config['dbname'])
    else:
        ## if database is home mysql database
        cnx = pymysql.connect(host=config_home.database_config['dbhost'],
                              port=config_home.database_config['port'],
                              user=config_home.database_config['dbuser'],
                              passwd=config_home.database_config['dbpass'],
                              db=config_home.database_config['dbname'])

    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    #sql_fathmm1 = "INSERT IGNORE impact (mutation_id, sample_id, tissue_type, cancer_type) VALUES(%s,%s,%s,%s)"
    sql_fathmm2 = "UPDATE impact SET fathmm_cancer_pred = %s, fathmm_cancer_score = %s WHERE mutation_id = %s;"

    count = 0
    entry = []
    f_list = []
    for k,v in results_dict.items():
        aa      = k
        prediction = v[0]
        score   = v[1]
        entry   = (prediction, score, aa)
        with cnx.cursor() as cursor:
            rows = cursor.execute(sql_fathmm2, entry)
            count += 1

    cnx.commit()
    cnx.close()

    return count





########## main ############

if __name__ == "__main__":

    results = fathmmResultsParser('fathmm_results.txt')
    #print(results)

    rows = fathmmInsert(results, 'home')
    print(rows)