#!/usr/bin python3

""" Cancer Genome Atlas Parser/Mutation table """

"""
Program:    gdc_parser
File:       gdc_parser.py
Version:    1.0
Date:       03.10.18
Function:   Parse gdc files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses gdc files for kinesin database mutation table insertion

Usage:
======
gdc_parser         SELF

Revision History:
=================
V1.2    22.01.19        PyMysql insert data, create dictionary      By: JJS
V1.3    23.01.19        Rewrite into function (insertMutation)          JJS     
V1.4    23.01.19        Added functions for source and impact           JJS
                        tables
"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import config_kinesin
import json

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

    sql_source = "INSERT INTO source_info (source_id, source_db, transcript_id, mutation_id) VALUES(%s,%s,%s,%s)"

    for x in mutations:
        if x[3] != "None": #  and x[4]=="missense_variant":
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
def parseMutations(file):
    """ Function to parse json files from Genomic Data Commons.
    Input                   file                file of JSON mutations
    Output                  mut_list  [0]       list of strings for mutation table
                            source_list  [1]    list of strings for source table
                            impact_list  [2]    list of strings for impact table
    """

    mut_list    = []
    mut_entry   = []
    source_entry= []
    source_list = []
    impact_entry = []
    impact_list = []

    source_db   = 'GDC'
    my_gene     = 'KIF11'

    ## first index is the item number, 'x'
    x = 0
    for k in file:
        try:
            gene_name       = str(mutations[x]['consequence'][0]['transcript']['gene']['symbol'])
            transcript      = str(mutations[x]['consequence'][0]['transcript']['gene']['gene_id'])
            protein         = str(mutations[x]['consequence'][0]['transcript']['aa_change'])
            consequence     = str(mutations[x]['consequence'][0]['transcript']['consequence_type'])
            genomic_id      = str(mutations[x]['genomic_dna_change'])
            ## eliminate 'g' in genomic_id
            genomic_id      = genomic_id.replace('g.', '')
            source_id       = str(mutations[x]['ssm_id'])
            mutation_type   = str(mutations[x]['mutation_subtype'])
            vep = str(mutations[x]['consequence'][0]['transcript']['annotation']['vep_impact'])
            if vep == '':
                vep = ' '
            polyphen = str(mutations[x]['consequence'][0]['transcript']['annotation']['polyphen_impact'])
            if polyphen == '':
                polyphen = ' '
            sift = str(mutations[x]['consequence'][0]['transcript']['annotation']['sift_impact'])
            if sift == '':
                sift = ' '
            x += 1
        except Error as e:
            print("Error", e)
        # put all strings into a list
        if gene_name == my_gene:            ## check that gene is KIF11
            mut_entry       = [genomic_id, 'y', ' ', mutation_type, consequence, protein, gene_name, 'homo sapiens']
            source_entry    = [source_id, source_db, transcript, protein]
            impact_entry = [protein, vep, sift, polyphen]
            impact_list.append(impact_entry)
            mut_list.append(mut_entry)
            source_list.append(source_entry)

    return mut_list, source_list, impact_list

# ******************************************************************************

########## main ############

if __name__ == "__main__":

    with open ('mutations.2018-10-03.json', 'r') as f:
        mutations = json.load(f)

    gdc_mut    = parseMutations(mutations)[0]
    gdc_source = parseMutations(mutations)[1]
    gdc_impact = parseMutations(mutations)[2]


    insertMutation(gdc_mut)
    insertSource(gdc_source)
    insertImpact(gdc_impact)

    i=0
    for k in gdc_mut:
        if k[5] != "None": # and k[4] == "missense_variant":
            i += 1
        #print(k)

    print('The total number of mutations is: ', i )

    f.close()


# how to find info for cds column?
# remember to delete all in mutation table on mysql database before running again

# checkout this tutorial: http://www.mysqltutorial.org/python-mysql-insert/
# make into function for inserting into table
# def main():
#     INSERT
#     INTO
#     mutation(genomic, coding, cds, mutation_type, consequence, protein, gene_name, organism)
#     					VALUES ('10:92616824C>T', 'Y', '1120', 'substitution', 'missense', 'L374F', 'KIF11', 'homo sapiens');
