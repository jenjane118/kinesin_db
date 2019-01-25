#!/usr/bin python3

""" Insertion module for kinesin database """

"""
Program:    dbinsert_module
File:       dbinsert_module.py
Version:    1.0
Date:       25.01.19
Function:   Insertion functions to populate kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program inserts attributes for tables in kinesin database.

Usage:
======
dbinsert_module         SELF

Revision History:
=================
V1.0    25.01.19        Initial version                             By: JJS

"""

# ******************************************************************************
# Import libraries

import sys
import pymysql
import config_kinesin
import mutation_parser

# ******************************************************************************
def insertMutation(mutations):
    """Inserts entry into mutation table.
    Input               mutations           list of mutations
    Output                                  inserted row into db
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_mutation = "INSERT INTO mutation (genomic, coding, cds, mutation_type, consequence, protein, " \
                   "gene_name, organism) VALUES(%s,%s,%s,%s,%s,%s,%s,%s)"

    for x in mutations:
        if x[5] != "None": #  and x[4]=="missense_variant":
            rows = cursor.execute(sql_mutation, x)

    cnx.commit()
    cnx.close()
# ******************************************************************************
def insertSource(mutations):
    """Inserts entry into source table.
    Input               mutations           list of mutations
    Output                                  inserted row into db
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_source = "INSERT INTO source_info (source_id, source_db, mutation_id) VALUES(%s,%s,%s)"

    for x in mutations:
        if x[2] != "None": #  and x[4]=="missense_variant":
            rows = cursor.execute(sql_source, x)

    cnx.commit()
    cnx.close()
# ******************************************************************************

def insertImpact(impact_list):
    """Inserts entry into impact table.
    Input               impact_list             list of entries for impact table
    Output                                      inserted row into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    # sql_impact = "INSERT INTO impact (mutation_id, custom_score, vep, sift_score, sift_prediction, " \
    #              "polyphen_score, polyphen_prediction, fathhm_score, fathhm_prediction, " \
    #              "clinvar_prediction) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    #
    sql_impact = "INSERT INTO impact (mutation_id, vep, sift_prediction, polyphen_prediction) " \
                 "VALUES(%s,%s,%s,%s)"


    for x in impact_list:
        if x[0] != "None":
            rows = cursor.execute(sql_impact, x)

    cnx.commit()
    cnx.close()

# ******************************************************************************

def insert_tissue(tissue_list):
    """Inserts entry into impact table.
    Input               tissue_list             list of entries for tissue table
    Output                                      inserted row into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_tissue = "INSERT INTO tissue (mutation_id, tissue_type, cancer_type) VALUES(%s,%s,%s)"

    for x in tissue_list:
        if x[0] != "None":
            rows = cursor.execute(sql_tissue, x)

    cnx.commit()
    cnx.close()

# ******************************************************************************
########## main ############

if __name__ == "__main__":

    gdc_att    = mutation_parser.parseGDC('KIF11', 'mutations.2018-10-03.json')
    gdc_mut    = gdc_att[0]
    gdc_source = gdc_att[1]
    gdc_impact = gdc_att[2]

    print(gdc_att[0])

    insertMutation(gdc_mut)
    insertSource(gdc_source)
    insertImpact(gdc_impact)