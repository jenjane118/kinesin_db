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
mutation text file and a function to parse the output file from VEP.

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
import db_query as q

#*****************************************************************************

def getVepInput(csv_file, results_file):

    veps = q.vepScores(gene)

    # run m.cosmicParser and extract gene_number and cds from dictionary
    line = ''
    with open(results_file, 'w') as outfile:
        mut_dict = m.cosmicParser(gene, csv_file)
        for mutation in veps:
            if mutation in mut_dict:
                number  = str(mut_dict[mutation][13])
                code    = str(mut_dict[mutation][3])
                line = number + ':' + code
                print(line, file=outfile)

    return outfile

# ******************************************************************************

def parseVep(my_gene, vep_file):
    ## parse vep results file for VEP impact, polyphen and sift scores/impact
    vep_entry       = []
    vep_impact_list = []
    # open csv file
    with open(vep_file, 'r') as file:
        vep_reader = csv.reader(file, delimiter='\t')
        for row in vep_reader:
            gene_name   = row[5]
            vep_impact  = row[4]
            sift        = row[27]
            polyphen    = row[28]
            # except Error as e:
            #     print("Error", e)
            if gene_name == my_gene:  ## check that gene is KIF11
                vep_entry = [vep_impact, sift, polyphen]
                vep_impact_list.append(vep_entry)

    return vep_impact_list
# *****************************************************************************

########## main ############

if __name__ == "__main__":

    gene = 'KIF11'
    cosmic_file = 'V87_38_MUTANT.csv'
    results_file = 'vep_input.txt'

    #getVepInput(cosmic_file, 'vep_input.txt')

    impact = parseVep(gene, 'vep_output.txt')

    for x in impact:
        print(x)


