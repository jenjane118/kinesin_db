#!/usr/bin python3

""" VEP access and Parser Module """

"""
Program:    vep_parse
File:       vep_parse.py
Version:    1.0
Date:       09.04.19
Function:   Parse VEP files for kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses VEP (Variant Effect Predictor) files (https://www.ensembl.org/Homo_sapiens/Tools/VEP) 
for impact results. It includes functions to create an input file for VEP webservice from the cosmic database
mutation text file and a function to parse the output file from VEP and update mysql database tables.

Usage:
======
vep_parse                SELF

Revision History:
=================
V1.0    09.04.19        Initial                                                 By: JJS
V1.2    10.06.19        Expand to parse more complete file using VEP website        JJS
                        Will start over to insert data in new impact table versus 
                        updating. 
v1.3    11.06.19        Finalised VEP parsing/insertion for all substitutions.      JJS
"""

#
#*****************************************************************************
# Import libraries

import csv
import mutation_parser as m
import re
import config_home
import pymysql
import config_kinesin

#*****************************************************************************

def vepScores(gene, database):
    """This function queries the kinesin database for list of mutations and formats
    for the VEP webservice. (https://www.ensembl.org/Tools/VEP)
    Input                   gene                       desired gene ('KIF11')
                            database                   kenobi or home database
    Output                  vep_list                   list of mutations in HGNS format
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

    mut_str = ''
    sub = ''
    nucl_sub = ''
    vep_list = []
    with cnx.cursor() as cursor:
        query = "SELECT genomic, cds FROM mutation WHERE mutation_type = 'substitution';"
        cursor.execute(query)
        temp = cursor.fetchall()

    for x in temp:
        genomic = str(x[0])
        genomic = genomic.replace('10:', '')
        sub = str(x[1])
        p = re.compile(r'(\w>\w)$')
        it = p.finditer(sub)
        for match in it:
            nucl_sub = match[0]

        mut_str = 'chr10:g.' + genomic + nucl_sub

        vep_list.append(mut_str)

    return vep_list

# *****************************************************************************

def parseVep2(my_gene, vep_file):
    """ Parse VEP flatfile results (downloaded after using webservice (https://www.ensembl.org/Multi/Tools/VEP?db=core)
    with file from getVepInput function) to get all predictions and scores for entry to kinesin
    database impact table.

    Input               my_gene                 Gene of interest (KIF11)
                        vep_file                flat file (results) downloaded from VEP website
    Output              impact_list             list of attributes for impact table update
    """

    vep2_entry       = []
    new_entry       = []
    vep2_impact_list = []
    # open csv file
    with open(vep_file, newline='') as file:
        vep_reader = csv.DictReader(file, delimiter='\t')
        for row in vep_reader:
            for k in row:
                if row[k] == '-':
                    row[k] = 'NA'
            try:
                gene_name       = row['SYMBOL']        #'KIF11'
                resnum          = row['Protein_position']
                aa_sub          = row['Amino_acids']
                if len(aa_sub)==1:
                    mutation = str(aa_sub + resnum + aa_sub)
                else:
                    q = re.compile(r'^(.*)/(.*)')
                    it = q.finditer(aa_sub)
                    for match in it:
                        aa1 = match.group(1)
                        aa2 = match.group(2)
                        mutation        = str(aa1+resnum+aa2)
                genomic         = row['#Uploaded_variation']        # mutation id (gene:genomic)
                # ## remove ENST gene info so can use to search db for cds attribute and pop missing info into impact table
                # p = re.compile(r'^chr')
                # genomic = p.sub('', genomic)
                # q = re.compile(r'g\.')
                # genomic = q.sub('', genomic)
                sift            = row['SIFT']
                if sift != 'NA':
                    s = re.compile(r'^(.*)\((.*)\)')
                    it = s.finditer(sift)
                    for match in it:
                        sift_pred = match.group(1)
                        sift_score = match.group(2)
                polyphen        = row['PolyPhen']
                if polyphen != 'NA':
                    it = s.finditer(polyphen)
                    for match in it:
                        polyphen_pred = match.group(1)
                        polyphen_score = match.group(2)
                condel          = row['Condel']
                if condel != 'NA':
                    it = s.finditer(condel)
                    for match in it:
                        condel = match.group(2)    # only interested in score
                cadd_raw        = row['CADD_raw']
                cadd_rank       = row['CADD_raw_rankscore']
                fathmm_rank     = row['FATHMM_converted_rankscore']
                fathmm_score    = row['FATHMM_score']
                fathmm_pred     = row['FATHMM_pred']
                metaSVM_pred    = row['MetaSVM_pred']
                metaSVM_rank    = row['MetaSVM_rankscore']
                metaSVM_score   = row['MetaSVM_score']
                mutpred_rank    = row['MutPred_rankscore']
                mutpred_score   = row['MutPred_score']
                mutassess_pred  = row['MutationAssessor_pred']
                mutassess_score = row['MutationAssessor_score']
                mutassess_rank  = row['MutationAssessor_score_rankscore']
                mutaster_rank   = row['MutationTaster_converted_rankscore']
                mutaster_pred   = row['MutationTaster_pred']
                mutaster_score  = row['MutationTaster_score']
                provean_rank    = row['PROVEAN_converted_rankscore']
                provean_score   = row['PROVEAN_score']
                provean_pred    = row['PROVEAN_pred']
                revel_rank      = row['REVEL_rankscore']
                revel_score     = row['REVEL_score']
            except NameError as e:
                print("Error", e)

            if gene_name == my_gene:  ## check that gene is KIF11
                vep2_entry = [mutation, sift_pred, sift_score, polyphen_pred, polyphen_score, condel,\
                              cadd_raw, cadd_rank,\
                              fathmm_score, fathmm_rank, fathmm_pred,\
                              metaSVM_score, metaSVM_rank, metaSVM_pred,\
                              mutpred_score, mutpred_rank,\
                              mutassess_score, mutassess_rank, mutassess_pred,\
                              mutaster_score, mutaster_rank, mutaster_pred,\
                              provean_score, provean_rank, provean_pred,\
                              revel_score, revel_rank]

                vep2_impact_list.append(vep2_entry)

    return vep2_impact_list

# *****************************************************************************

def updateImpact2(impact_list, database):
    """ Connects to kinesin database. Updates relevant entries in impact table with all metrics from VEP.
    First must find mutation_id for each entry based on cds attribute from impact_list. Use this to update impact table.
    Input                   impact_list                     List of attributes for update of impact table
                            database                        home or kenobi
    Output                  i                               number of updated rows
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

    sql_impact = "INSERT IGNORE impact (mutation_id, polyphen_pred, polyphen_score, sift_pred, sift_score,"\
                    "condel, cadd_score, cadd_rank, fathmm_score, fathmm_rank, fathmm_pred,"\
                    "metaSVM_score, metaSVM_rank, metaSVM_pred, mutpred_score, mutpred_rank,"\
                    "mutassessor_score, mutassessor_rank, mutassessor_pred, muttaster_score, muttaster_rank, muttaster_pred,"\
                    "provean_score, provean_rank, provean_pred, revel_score, revel_rank)\
                     VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    i = 0
    for x in impact_list:
        with cnx.cursor() as cursor:
            rows = cursor.execute(sql_impact, x)
            i += 1
    cnx.commit()
    cnx.close()

    return i

# ******************************************************************************

########## main ############

if __name__ == "__main__":

    gene = 'KIF11'

    #input_list = vepScores(gene, 'home')

    #with open ('vep_input.txt', 'w') as outfile:
    #    for x in input_list:
    #       print(x, file=outfile)

    # use mutations listed in file in vep webservice (or use command line request)
    # output from webservice in 'vep_complete_results.txt'

    impact = parseVep2(gene, 'vep_complete_results.txt')
    number = updateImpact2(impact, 'home')
    print(number)


