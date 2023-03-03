#!/usr/bin python3

""" Domain mapping function """

"""
Program:    domain_mapper
File:       domain_mapper.py
Version:    1.0
Date:       29.01.19
Function:   Maps amino acid residue to appropriate kinesin domain
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program includes a function to map amino acid residues onto pfam kinesin domain architecture
https://pfam.xfam.org/protein/P52732.

Usage:
======
domain_mapper           SELF

Revision History:
=================
V1.0    29.01.19        Initial version                             By: JJS
V1.1    13.07.19        Change domain names                             JJS
V1.2    17.08.19        Made separate function for IDRs                 JJS
"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import config_home

# ******************************************************************************

def getMutations(database):

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

    mutation_list = []

    # query list of mutations:  aa_num, sample_id, tissue_type
    query = "SELECT protein FROM mutation m  WHERE m.consequence = 'missense';"
            #"AND t.tissue_type in ('Rectum', 'Colon', 'large_intestine')"
            #"AND m.protein = t.mutation_id;"

    with cnx.cursor() as cursor:
        cursor.execute(query)
        temp = cursor.fetchall()

    # print string for each mutation to text file
    with open ('missense.txt', 'w') as outfile:
        mutation_string = ''

        for x in temp:
            #residue_list = temp.read().splitlines()
            #for row in temp:

            mutation_string = str(x[0])
            mutation_list.append(mutation_string)
        for y in mutation_list:
            print(y, end=' ', file=outfile)


    return mutation_list




# ******************************************************************************

def domainMapper(mutation):
    """This function maps the amino acid mutation to a domain on the KIF11-Human pfam architecture
    (Using uniprot id P52732).

    Input               mutation            mutated amino acid residue
    Output              domain_name         name of relevant domain
                                            
    """
    import re

    domain_dict = {}
    for x in range(1,2000):
        if x in range(24,359):
            domain_dict[str(x)] = 'kinesin motor'
        elif x in range(916, 1053): #> 915 and x < 1054:
            domain_dict[str(x)] = 'tail-binding'
        else:
            domain_dict[str(x)] = 'coiled-coil/disorder'

    #find residue number
    residue_num = ''
    r = re.compile(r'[A-Z]+(\d+)[^0-9]')
    it = r.finditer(mutation)
    for match in it:
        residue_num = str(match.group(1))

    if residue_num in domain_dict:
        return  domain_dict[residue_num]
    else:
        return  'UNK'

# *************************************************************************************

def domainMapper2(mutation):
    """This function maps the amino acid mutation to a domain on the KIF11-Human pfam architecture
    (using uniprot id P52732) and includes idr's (intrinsically disordered regions can include
    important functional sites such as phosphorylation sites or ppi mediators/regulators. Identified
    with e-Driver suppl data (used Foldindex)(http://github.com/eduardporta/e-Driver.git).

    Input               mutation            mutated amino acid residue
    Output              domain_name         name of relevant domain

    """
    import re

    domain_dict = {}
    for x in range(1, 2000):
        if x in range(24, 359):
            domain_dict[str(x)] = 'kinesin motor'
        elif x in range(916, 1053):  # > 915 and x < 1054:
            domain_dict[str(x)] = 'tail-binding'
        elif x in range(751, 885):
            domain_dict[str(x)] = 'idr1'
        elif x in range(422, 496):
            domain_dict[str(x)] = 'idr2'
        else:
            domain_dict[str(x)] = 'coiled-coil/disorder'

    # find residue number
    residue_num = ''
    r = re.compile(r'[A-Z]+(\d+)[^0-9]')
    it = r.finditer(mutation)
    for match in it:
        residue_num = str(match.group(1))

    if residue_num in domain_dict:
        return domain_dict[residue_num]
    else:
        return 'UNK'


# *************************************************************************************

########## main ############

if __name__ == "__main__":

    with open('missense.txt', 'r') as file:
        residue_list = file.read().splitlines()
        print(residue_list)
        idr1_list = []
        idr2_list = []
        for row in residue_list:
            p = domainMapper2(row)
            if p == 'idr1':
                idr1_list.append(row)
            elif p == 'idr2':
                idr2_list.append(row)
            else:
                pass
        print(len(idr1_list))
        print(len(idr2_list))
    with open('idr1_mutations.txt', 'w') as outfile:
        for x in idr1_list:
            print(str(x), end=' ', file=outfile)
