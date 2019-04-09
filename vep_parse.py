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
import re

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
## match up mutation with impact info




########## main ############

if __name__ == "__main__":

    gene = 'KIF11'
    cosmic_file = 'V87_38_MUTANT.csv'
    results_file = 'vep_input.txt'

    #getVepInput(cosmic_file, 'vep_input.txt')

    impact = parseVep(gene, 'vep_output.txt')





    ## zip together vep_input file and impact list
    impact_list = []

    with open ('vep_input.txt', 'r') as file:
        for line in file:
            gene_list = file.read().splitlines()
    impact_list = list(zip(gene_list, impact))

    #reformat to be a 4 attribute list for each entry
    new_list = []
    new_entry = []
    sift_pred = ''
    sift_score = ''
    for entry in impact_list:
        mutation    = str(entry[0])
        ## remove ENST gene info so can use to search db and pop in missing info
        p           = re.compile(r'^ENST00000260731:c\.')
        mutation    = p.sub('', mutation)
        vep         = str(entry[1][0])

        ## need to parse out score and prediction
        polyphen    = str(entry[1][1])
        if polyphen == '-':
            poly_pred = 'UNK'
            poly_score = 'UNK'
        else:
            q       = re.compile(r'^(.*)\((.*)\)')
            it      = q.finditer(polyphen)
            for match in it:
                poly_pred  = match.group(1)
                poly_score = match.group(2)

        sift        = str(entry[1][2])
        if sift     == '-':
            sift_pred  = 'UNK'
            sift_score = 'UNK'
        else:
            s      = re.compile(r'^(.*)\((.*)\)')
            it      = s.finditer(sift)
            for match in it:
                sift_pred   = match.group(1)
                sift_score  = match.group(2)

        new_entry   = [mutation, vep, poly_pred, poly_score, sift_pred, sift_score]
        new_list.append(new_entry)
    print(new_list)