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
def insertCosmicSource(mutation_dict):
    """Inserts entries from Cosmic mutations into source_info table.
    Input               mutation_dict         dictionary of attributes from Cosmic csv file
    Output                                  inserted row into db source_info table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_source = "INSERT INTO source_info (source_id, source_db, mutation_id) VALUES(%s,%s,%s)"

    source_line = []
    source_list = []
    # extract values from mutation_dict and put in list
    for k in mutation_dict:
        if k != "None":
            source_line = [mutation_dict[k][8], mutation_dict[k][7], k]
            source_list.append(source_line)

    # use list to populate source table
    for x in source_list:
            rows = cursor.execute(sql_source, x)

    cnx.commit()
    cnx.close()

# ******************************************************************************
def insertMutation(mutations):
    """Inserts entry into mutation table.
    Input               mutations           list of mutations from gdc
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
                   "gene_name, organism, domain) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    for x in mutations:
        if x[5] != "None": #  and x[4]=="missense_variant":
            rows = cursor.execute(sql_mutation, x)

    cnx.commit()
    cnx.close()
# ******************************************************************************
def insertSource(mutations):
    """Inserts entry into source table.
    Input               mutations           list of mutations from gdc
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

def insertImpact(impact):
    """Inserts entry into impact table.
    Input               impact                  list of entries for impact table
    Output                                      inserted row into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_impact = "INSERT INTO impact (mutation_id, vep, sift_prediction, polyphen_prediction, " \
                 "fathhm_score, fathhm_prediction) VALUES(%s,%s,%s,%s,%s,%s)"
    i = 0
    for x in impact:
        if x[0] != "None":
            rows = cursor.execute(sql_impact, x)
            i += 1
    cnx.commit()
    cnx.close()
    return  i
# ******************************************************************************

def insertTissue(mutation_dict):
    """Inserts entry into impact table.
    Input               mutation_dict               dict of mutation attributes from cosmic
    Output                                          inserted row into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_tissue = "INSERT INTO tissue (mutation_id, tissue_type, cancer_type) VALUES(%s,%s,%s)"

    tissue_line = []
    tissue_list = []
    # extract values from mutation_dict and put in list
    for k in mutation_dict:
        if k != "None":
            tissue_line = [k, mutation_dict[k][11], mutation_dict[k][12]]
            tissue_list.append(tissue_line)

    # use list to populate tissue table
    for x in tissue_list:
        rows = cursor.execute(sql_tissue, x)

    cnx.commit()
    cnx.close()

# ******************************************************************************

def insertCosmicMutation(mutation_dict):
    """Inserts entry into mutation table. (have to do each table in sep function to avoid foreign key constraints)
    Input               mutation_dict       dictionary of attributes from Cosmic csv file
    Output                                  inserts attributes into mutation table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    # insert mutation attributes from cosmic, will replace earlier entry with more complete info from cosmic
    sql_mutation    = "REPLACE INTO mutation (genomic, coding, cds, mutation_type, consequence, protein, " \
                         "gene_name, organism, domain) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    mutation_entry  = []
    mutation_list   = []

    for k in mutation_dict:
        if k != "None":
            mutation_entry = [mutation_dict[k][0], mutation_dict[k][1], mutation_dict[k][2], mutation_dict[k][3],
                              mutation_dict[k][4], k, mutation_dict[k][13], mutation_dict[k][5],
                              mutation_dict[k][6]]
            # insert cosmic impact info for mutations not in gdc
            mutation_list.append(mutation_entry)

    i = 0
    # use list to populate mutation table
    for x in mutation_list:
        rows = cursor.execute(sql_mutation, x)
        i += 1

    cnx.commit()
    cnx.close()

    return i

# ******************************************************************************
def insertCosmicImpact(mutation_dict):
    """Inserts entry into mutation table. (have to do each table in sep function to avoid foreign key constraints)
    Input               mutation_dict       dictionary of attributes from Cosmic csv file
    Output                                  inserts attributes into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    # insert impact attributes for mutations in cosmic but ignore mutations already in db
    sql_impact      = "INSERT IGNORE impact (mutation_id, fathhm_score, fathhm_prediction) VALUES(%s,%s,%s)"

    impact_entry    = []
    impact_list     = []

    for k in mutation_dict:
        if k != "None":
            impact_entry   = [k, mutation_dict[k][9], mutation_dict[k][10]]
            impact_list.append(impact_entry)
    i = 0
    for y in impact_list:
        rows = cursor.execute(sql_impact, y)
        i += 1

    cnx.commit()
    cnx.close()

    return i

# ******************************************************************************
def insertCosmicTissue(mutation_dict):
    """Inserts entry into mutation table. (have to do each table in sep function to avoid foreign key constraints)
    Input               mutation_dict       dictionary of attributes from Cosmic csv file
    Output                                  inserts attributes into tissue table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_tissue      = "INSERT INTO tissue (mutation_id, tissue_type, cancer_type) VALUES(%s,%s,%s)"

    tissue_line     = []
    tissue_list     = []

    for k in mutation_dict:
        if k != "None":
            tissue_line = [k, mutation_dict[k][11], mutation_dict[k][12]]
            tissue_list.append(tissue_line)

    i = 0

    # use list to populate tissue table
    for z in tissue_list:
        rows = cursor.execute(sql_tissue, z)
        i += 1

    cnx.commit()
    cnx.close()

    return i
# ******************************************************************************

########## main ############

if __name__ == "__main__":

    gdc_att    = mutation_parser.parseGDC('KIF11', 'mutations.2018-10-03.json')
    gdc_mut    = gdc_att[0]
    gdc_source = gdc_att[1]
    gdc_impact = gdc_att[2]

    cosmic_mutation   = mutation_parser.cosmicParser('KIF11', 'V87_38_MUTANT.csv')
    impact_all = mutation_parser.combineImpact('KIF11', 'mutations.2018-10-03.json', 'V87_38_MUTANT.csv')
    #print(impact_all)
    insertMutation(gdc_mut)
    insertCosmicMutation(cosmic_mutation)
    insertSource(gdc_source)
    insertCosmicSource(cosmic_mutation)
    insertImpact(impact_all)
    #insertImpact(gdc_impact) (need a seperate function)
    insertCosmicImpact(cosmic_mutation)
    insertCosmicTissue(cosmic_mutation)

