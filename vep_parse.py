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
This program parses VEP (Variant Effect Predictor) files (https://www.ensembl.org/info/docs/tools/vep/index.html) 
for missing impact results. It includes functions to create an input file for VEP webservice from the cosmic database
mutation text file and a function to parse the output file from VEP and update mysql database tables.

Usage:
======
vep_parse                SELF

Revision History:
=================
V1.0    09.04.19        Initial                                     By: JJS

"""

#
#*****************************************************************************
# Import libraries

import csv
import mutation_parser as m
import re
import config_home
import pymysql

#*****************************************************************************

def vepScores(gene):
    """This function queries the kinesin database for list of mutations that need vep scores.
    Input                   gene                       desired gene ('KIF11')
    Output                  vep_list
    """

    # Connect to MySQL Database
    cnx = pymysql.connect(host=config_home.database_config['dbhost'],
                          port=config_home.database_config['port'],
                          user=config_home.database_config['dbuser'],
                          passwd=config_home.database_config['dbpass'],
                          db=config_home.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    vep_dict = {}
    vep_list = []
    with cnx.cursor() as cursor:
        query = "SELECT mutation_id FROM impact WHERE vep = 'UNK';"
        cursor.execute(query)
        temp = cursor.fetchall()

    for x in temp:
        vep_list.append(x[0])

    return vep_list

#*****************************************************************************

def getVepInput(csv_file, results_file):
    """ This function runs mutation_parser.cosmicParser function to extract 'gene number' and cds from
    cosmic download file to make a list of mutations in correct format to submit to VEP webservices
    (https://www.ensembl.org/Multi/Tools/VEP?db=core).

    Input                   csv_file                flat file COSMIC mutations (downloaded)
                            results_file            name of outfile for submitting to VEP
    Output                  outfile                 file of mutations in correct format for submission to webservices
    """

    # run vepScores function to query kinesin database for entries lacking VEP predictions
    veps = vepScores(gene)

    # run m.cosmicParser and extract gene_number and cds from resulting mutation dictionary to make list of mutations
    line = ''
    with open(results_file, 'w') as outfile:
        mut_dict = m.cosmicParser(gene, csv_file)
        for mutation in veps:
            if mutation in mut_dict:
                number  = str(mut_dict[mutation][13])       #'gene number' = ensembl id (ENST00000260731)
                code    = str(mut_dict[mutation][3])        # cds
                line = number + ':' + code
                print(line, file=outfile)

    return outfile

# ******************************************************************************

def parseVep(my_gene, vep_file):
    """ Parse VEP flatfile results (downloaded after using webservice (https://www.ensembl.org/Multi/Tools/VEP?db=core)
    with file from getVepInput function) to get VEP, sift and polyphen predictions and scores for entry to kinesin
    database impact table.

    Input               my_gene                 Gene of interest (KIF11)
                        vep_input               file of mutations submitted to VEP
                        vep_file                flat file (results) downloaded from VEP website
    Output              impact_list             list of attributes for impact table update
    """

    vep_entry       = []
    vep_impact_list = []
    # open csv file
    with open(vep_file, 'r') as file:
        vep_reader = csv.reader(file, delimiter='\t')
        for row in vep_reader:
            gene_name   = row[5]        #'KIF11'
            mutation    = row[0]        # mutation id (gene:cds)
            vep_impact  = row[4]        # prediction only
            sift        = row[27]       # includes prediction and score in one entry
            polyphen    = row[28]       # incl prediction and score in one entry
            # except Error as e:
            #     print("Error", e)
            if gene_name == my_gene:  ## check that gene is KIF11
                vep_entry = [mutation, vep_impact, sift, polyphen]
                vep_impact_list.append(vep_entry)

    # reformat to be a 4 attribute list for each entry suitable for updating impact table
        impact_list = []
        new_entry = []

        for entry in vep_impact_list:
            mutation = str(entry[0])
            ## remove ENST gene info so can use to search db for cds attribute and pop missing info into impact table
            p = re.compile(r'^ENST00000260731:c\.')
            mutation = p.sub('', mutation)
            vep = str(entry[1])
            ## need to parse out score and prediction from combined entry in vep results file
            polyphen = str(entry[2])
            if polyphen == '-':
                poly_pred = 'UNK'
                poly_score = 'UNK'
            else:
                q = re.compile(r'^(.*)\((.*)\)')
                it = q.finditer(polyphen)
                for match in it:
                    poly_pred = match.group(1)
                    poly_score = match.group(2)

            sift = str(entry[3])
            if sift == '-':
                sift_pred = 'UNK'
                sift_score = 'UNK'
            else:
                s = re.compile(r'^(.*)\((.*)\)')
                it = s.finditer(sift)
                for match in it:
                    sift_pred = match.group(1)
                    sift_score = match.group(2)

            new_entry = [mutation, vep, poly_pred, poly_score, sift_pred, sift_score]
            impact_list.append(new_entry)

    return impact_list

# *****************************************************************************

def updateImpact(impact_list):
    """ Connects to kinesin database. Updates relevant entries in impact table with VEP, polyphen and sift predictions.
    First must find mutation_id for each entry based on cds attribute from impact_list. Use this to update impact table.
    Input                   impact_list                     List of attributes for update of impact table
    Output                  i                               number of updated rows
    """

    # Connect to MySQL Database
    cnx = pymysql.connect(host=config_home.database_config['dbhost'],
                          port=config_home.database_config['port'],
                          user=config_home.database_config['dbuser'],
                          passwd=config_home.database_config['dbpass'],
                          db=config_home.database_config['dbname'])
    cursor = cnx.cursor(pymysql.cursors.DictCursor)

    insert_list = []

    # use cds to obtain mutation_id from mutation table
    query = "SELECT protein FROM mutation WHERE cds = %s ;"
    for x in impact_list:
        mutation = str.format(x[0])         # retrieve mutation (cds) from parsed impact list
        with cnx.cursor() as cursor:
            cursor.execute(query, (mutation))
            temp = cursor.fetchone()

        # inserts amino acid change as first element of list
        x.insert(6, temp[0])
        # remove cds from end of list
        x.pop(0)
        insert_list.append(x)
    #print(insert_list)

    sql_impact = "UPDATE impact SET vep = %s , polyphen_prediction = %s , polyphen_score =  %s , sift_prediction =  %s ,"\
                    "sift_score =  %s  WHERE mutation_id = %s;"
    i = 0
    for x in insert_list:
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
    cosmic_file = 'V87_38_MUTANT.csv'


    input_file = getVepInput(cosmic_file, 'vep_input2.txt')

    impact = parseVep(gene, 'VEP_output2.txt')
    print(impact)

    number = updateImpact(impact)
    print(number)
