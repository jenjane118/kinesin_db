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
This program calls on parsing module functions and inserts attributes for tables in
kinesin database.

Usage:
======
dbinsert_module         SELF

Revision History:
=================
V1.0    25.01.19        Initial version                             By: JJS
V1.1    26.01.19        Complete insertion from both databases          JJS

"""

# ******************************************************************************
# Import libraries

import sys
import pymysql
import config_kinesin
#import config_home
import mutation_parser


# ******************************************************************************
def insertCosmicSource(mutation_dict):
    """Inserts entries from Cosmic mutations into source_info table.
    Input               mutation_dict         dictionary of attributes from Cosmic csv file
    Output                                  inserted row into db source_info table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host  =config_kinesin.database_config['dbhost'],
                          #port  =config_kinesin.database_config['port'],
                          user  =config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db    =config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_source = "INSERT INTO source_info (source_id, source_db, mutation_id) VALUES(%s,%s,%s)"

    source_line = []
    source_list = []
    # extract values from mutation_dict and put in list
    for k in mutation_dict:
        if k != "None":
            source_line = [mutation_dict[k][9], mutation_dict[k][8], k]
            source_list.append(source_line)

    # use list to populate source table
    for x in source_list:
            rows = cursor.execute(sql_source, x)

    cnx.commit()
    cnx.close()

# ******************************************************************************
def insertMutation(mutations):
    """Inserts entry into mutation table.
    Input               mutations           list of attributes from gdc
    Output                                  inserted row into db
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_mutation = "INSERT INTO mutation (protein, resnum, genomic, coding, cds, mutation_type, consequence, " \
                   "gene_name, organism, domain) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    gdc_mut = mutations[0]
    for x in gdc_mut:
        if x[0] != "None": #  and x[4]=="missense_variant":
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
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_source = "INSERT INTO source_info (source_id, source_db, mutation_id) VALUES(%s,%s,%s)"
    gdc_source = mutations[1]
    for x in gdc_source:
        if x[2] != "None": #  and x[4]=="missense_variant":
            rows = cursor.execute(sql_source, x)

    cnx.commit()
    cnx.close()
# ******************************************************************************

def insertCombinedImpact(impact):
    """Inserts entry into impact table.
    Input               impact                  list of combined entries for impact table
    Output                                      inserted row into impact table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_impact = "INSERT INTO impact (mutation_id, vep, sift_prediction, polyphen_prediction, " \
                 "fathmm_score, fathmm_prediction) VALUES(%s,%s,%s,%s,%s,%s)"
    i = 0
    for x in impact:
        if x[0] != "None":
            rows = cursor.execute(sql_impact, x)
            i += 1
    cnx.commit()
    cnx.close()
    return  i

# ******************************************************************************

def insertGdcImpact(mutations):
    """Inserts entry into impact table.
    Input               impact                  list of mutation attributes from gdc
    Output                                      inserted row into impact table
                        i                       number of inserted rows
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_impact = "INSERT IGNORE impact (mutation_id, vep, sift_prediction, polyphen_prediction) VALUES(%s,%s,%s,%s)"
    gdc_impact = mutations[2]
    i = 0
    for x in gdc_impact:
        if x[0] != "None":
            rows = cursor.execute(sql_impact, x)
            i += 1
    cnx.commit()
    cnx.close()
    return  i

# ******************************************************************************

def insertCosmicTissue(cos_tissue_list):
    """Inserts entry into tissue table.
    Input               mutation_tissue             list of mutation attributes from cosmic
    Output                                          inserted row into tissue table
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_tissue = "INSERT IGNORE tissue (mutation_id, sample_id, tissue_type, cancer_type) VALUES(%s,%s,%s,%s)"

    tissue_line = []
    tissue_list = []
    # use list to populate tissue table
    for x in cos_tissue_list:
        if x[0] != "None":
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
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    # insert mutation attributes from cosmic, will replace earlier entry with more complete info from cosmic
    sql_mutation    = "REPLACE INTO mutation (protein, resnum, genomic, coding, cds, mutation_type, consequence, " \
                         "gene_name, organism, domain) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    mutation_entry  = []
    mutation_list   = []

    for k in mutation_dict:
        if k != "None":
            mutation_entry = [k, mutation_dict[k][0], mutation_dict[k][1], mutation_dict[k][2], mutation_dict[k][3],
                              mutation_dict[k][4], mutation_dict[k][5], mutation_dict[k][12],
                              mutation_dict[k][6], mutation_dict[k][7]]
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
                          #port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    # insert impact attributes for mutations in cosmic but ignore mutations already in db
    sql_impact      = "INSERT IGNORE impact (mutation_id, fathmm_score, fathmm_prediction) VALUES(%s,%s,%s)"

    impact_entry    = []
    impact_list     = []

    for k in mutation_dict:
        if k != "None":
            impact_entry   = [k, mutation_dict[k][11], mutation_dict[k][10]]
            impact_list.append(impact_entry)
    i = 0
    for y in impact_list:
        rows = cursor.execute(sql_impact, y)
        i += 1

    cnx.commit()
    cnx.close()

    return i


# ******************************************************************************

def insertGdcTissue(tissue_list):
    """Inserts entry into impact table.
    Input               tissue_list             list of mutation attributes from gdc for tissue
    Output                                      inserted row into impact table
                        i                       number of inserted rows
    """
    # Connect to MySQL database
    cnx = pymysql.connect(host=config_kinesin.database_config['dbhost'],
                          port=config_kinesin.database_config['port'],
                          user=config_kinesin.database_config['dbuser'],
                          passwd=config_kinesin.database_config['dbpass'],
                          db=config_kinesin.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    sql_tissue = "INSERT IGNORE tissue (mutation_id, sample_id, tissue_type, cancer_type) VALUES(%s,%s,%s,%s)"

    i = 0
    for x in tissue_list:
        if x[0] != "None":
            rows = cursor.execute(sql_tissue, x)
            i += 1
    cnx.commit()
    cnx.close()
    return  i

# ******************************************************************************

########## main ############

if __name__ == "__main__":

    gdc_att         = mutation_parser.parseGDC('KIF11', 'mutations.2019-01-23.json')
    cosmic_mutation = mutation_parser.cosmicParser('KIF11', 'V87_38_MUTANT.csv')
    impact_all      = mutation_parser.combineImpact('KIF11', 'mutations.2019-01-23.json', 'V87_38_MUTANT.csv')
    cos_tissue_list = mutation_parser.tissueCosmic('KIF11', 'V87_38_MUTANT.csv')
    tissueGDC       = mutation_parser.tissueGDC('KIF11', 'results.json')

    # insert commands must be in this order
    insertMutation(gdc_att)
    insertCosmicMutation(cosmic_mutation)
    insertSource(gdc_att)
    insertCosmicSource(cosmic_mutation)
    insertCombinedImpact(impact_all)
    insertGdcImpact(gdc_att)
    insertCosmicImpact(cosmic_mutation)

    insertCosmicTissue(cos_tissue_list)

    insertGdcTissue(tissueGDC)