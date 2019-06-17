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
import config_kinesin
import Bio.Data.IUPACData

#*****************************************************************************

def parseClinvar(gene, csv_file):
    """This function parses Clinvar results file.
    Input               gene                KIF11
                        csv_file            file of Clinvar results returned from web query
    Output              cv_list             list of mutation identifier and clinical significance
    """
    #new_mutation = ''
    #mut = ''
    #clin_sig = ''
    cv_entry = []
    cv_list  = []
    # open csv file
    with open(csv_file, 'r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        #first_row = True
        for row in csv_reader:
            clin_sig = row['Clinical significance (Last reviewed)']  # clinical significance
            # remove review date in parentheses
            r = re.compile(r'\(.*\)$')
            clin_sig = r.sub('', clin_sig)
            mut_name    = row['Name']
                    ##NM_004523.3(KIF11):c.2T>C (p.Met1Thr)
            p       = re.compile(r'.+(\(KIF11\)).+p.(\w+)')
            it      = p.finditer(mut_name)
            for match in it:
                mut = str(match.group(2))
                # extract 3-letter amino acid codes and resnum
                q       = re.compile(r'(\D+)(\d+)(\D{3})')
                ob      = q.finditer(mut)
                for x in ob:
                    old_aa3 = str(x.group(1))
                    resnum = str(x.group(2))
                    new_aa3 = str(x.group(3))
                # reformat protein entry for one letter IUPAC code
                    try:
                        old_aa1 = Bio.Data.IUPACData.protein_letters_3to1[old_aa3]
                    except KeyError:
                        old_aa1 = '*'
                    try:
                        new_aa1 = Bio.Data.IUPACData.protein_letters_3to1[new_aa3]
                    except KeyError:
                        new_aa1 = '*'
                    new_mutation = old_aa1 + resnum + new_aa1
            #cv_list.append(new_mutation)
                    cv_entry = [clin_sig, new_mutation]
                    cv_list.append(cv_entry)
    return cv_list

# ******************************************************************************

def clinvarUpdate(clinvar_list, database):
    """ This function uses list of clinvar attributes to update impact table in kinesin database.
    Input               clinvar_list                    list of attributes (genomic location, clinical significance)
                        database                        which database used
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

    sql_clinvar = "UPDATE impact SET clinvar_prediction = %s WHERE mutation_id = %s;"

    count = 0
    for item in clinvar_list:
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

    db = 'home'

    results = parseClinvar(mygene, myfile)
    print(results)

    successful = clinvarUpdate(results, db)
    print(successful)