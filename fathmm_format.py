#!/usr/bin python3

""" FATHMM Formatting Script"""

"""
Program:    fathmm_format
File:       fathmm_format.py
Version:    1.0
Date:       13.07.19
Function:   Format mutations for use with FATHMM webserver
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
Format mutations for use in FATHMM webserver at http://fathmm.biocompute.org.uk/cancer.html

Usage:
======
fathmm_format                SELF

Revision History:
=================
V1.0    13.07.19        Initial                                     By: JJS

"""

#
#*****************************************************************************
# Import libraries

import config_home
import pymysql
import config_kinesin

#*****************************************************************************

def fathmm_format(gene, database):
    """This function queries the kinesin database for list of mutations and formats
    for the FATHMM webservice (http://fathmm.biocompute.org.uk/cancer.html) by adding
    commas to list of mutations.

    Input                   gene                       desired gene ('KIF11')
                            database                   kenobi or home database
    Output                  fathmm_list                list of mutations in simple comma separated format
                                                       (e.g. 'P52732   N342V, E101V')
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

    fathmm_str = ''
    fathmm_list = []

    with cnx.cursor() as cursor:
        query = "SELECT protein FROM mutation WHERE mutation_type = 'substitution';"
        cursor.execute(query)
        temp = cursor.fetchall()

    for x in temp:
        fathmm_list.append(str(x[0]))
    separator = ','
    fathmm_str = separator.join(fathmm_list)
    fathmm_str = 'P52732 ' + fathmm_str
    with open ('fathmm_mutations.txt', 'w') as outfile:
        print(fathmm_str, file=outfile)






    return outfile


# *****************************************************************************
########## main ############

if __name__ == "__main__":

    fathmm_format('KIF11', 'home')
