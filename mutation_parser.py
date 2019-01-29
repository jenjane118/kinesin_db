#!/usr/bin python3

""" Cancer Genome Atlas Parser/Mutation table """

"""
Program:    mutation_parser
File:       mutation_parser.py
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
This program parses gdc files (https://portal.gdc.cancer.gov) for kinesin database population.

Usage:
======
mutation_parser         SELF

Revision History:
=================
V1.2    22.01.19        PyMysql insert data, create dictionary      By: JJS
V1.3    23.01.19        Rewrite into functions                          JJS     
V1.4    23.01.19        Added functions for source and impact           JJS
                        tables
V1.5    25.01.19        Added file opening and gene specification to    JJS
                        parseMutation function. Removed insert          
                        functions to a new module. Renamed module 
                        (previously gdc_parser).
"""

# ******************************************************************************
# Import libraries

import re
import sys
import pymysql
import config_kinesin
import json
import csv


# ******************************************************************************
def parseGDC(gene, json_file):
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

    with open (json_file) as f:
        mutations = json.load(f)
        ## first index is the item number, 'x'
        x = 0
        for k in mutations:
            try:
                gene_name       = str(mutations[x]['consequence'][0]['transcript']['gene']['symbol'])
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
                coding          = 'y'
                cds             = ' '
                organism        = 'Homo sapiens'
                domain          = ' '              # determine later
                x += 1
            except Error as e:
                print("Error", e)

            # make terminology more consistent between databases
            p = re.compile(r'_variant$')
            consequence = p.sub('', consequence)
            q = re.compile(r'^Single base')
            mutation_type = q.sub('', mutation_type)

            # put all strings into lists for each table
            if gene_name == my_gene:            ## check that gene is KIF11
                mut_entry       = [genomic_id, coding, cds, mutation_type, consequence, protein, gene_name, organism, domain]
                source_entry    = [source_id, source_db, protein]
                impact_entry    = [protein, vep, sift, polyphen]
                impact_list.append(impact_entry)
                mut_list.append(mut_entry)
                source_list.append(source_entry)
    f.close()

    return mut_list, source_list, impact_list

# ******************************************************************************
def cosmicParser(my_gene, csv_file):
    """Function to parse mutation files from csv files from COSMIC database
    Input               file                    csv file download
    Output              id_list                 list of mutation_id's
                        mutation_dict           dictionary of attributes of mutation and impact

    """

    mut_entry       = []
    mutation_dict   = {}

    # open csv file
    with open(csv_file) as file:
        csv_reader = csv.reader(file, delimiter=',')
        # parse out info into attribute names
        for row in csv_reader:
            try:
                gene_name       = row[0]
                genomic_id      = row[1]
                coding          = 'y'
                protein         = row[18]
                # eliminate 'p.'
                protein         = protein.replace('p.', '')
                cds             = row[17]
                description     = row[19]           # parse out type and consequence below
                organism        = 'Homo sapiens'
                domain          = ' '               # calculate later

                # source_info table
                source_db       = 'COSMIC'
                source_id       = row[16]
                # impact table
                fathhm_score    = row[28]
                fathhm_pred     = row[27]

                # tissue table
                tissue_type     = row[7]
                cancer_type     = row[12]
            except Error as e:
                print("Error", e)

            # use regex to parse mutation type and consequence from overall mutation description
            p  = re.compile(r'(^\w+)\s-\s(.+$)')
            it = p.finditer(description)
            for match in it:
                mutation_type = match.group(1)  # need all before hyphen
                consequence   = match.group(2)  # after the hyphen
                # makes more consistent terminology between databases
                consequence = consequence.replace('coding silent', 'synonymous')

            # make dictionary to fill in missing attribute fields in mutation table before insertion
            if protein != '?' and gene_name == my_gene:  ## check that gene is KIF11, don't include ? for aa_change
                mutation_dict[protein] = (genomic_id, coding, cds, mutation_type, consequence, organism, domain,
                                          source_db, source_id, fathhm_score, fathhm_pred, tissue_type,
                                          cancer_type, gene_name)
    file.close()
    return mutation_dict

# *************************************************************************************
def combineImpact(gene, json_file, csv_file):
    """Creates list of combined impact attributes with attributes from GDC .json and Cosmic .csv files.
    Input           gene                gene name
                    json_file           gdc file in json format
                    csv_file            cosmic file in csv format
    Output          total_impact        updated impact list for insertion
    """
    gdc_att = parseGDC(gene, json_file)
    cosmic_mutation = cosmicParser(gene, csv_file)
    gdc_impact = gdc_att[2]
    # add fathhm attributes from cosmic file to impact list
    y = []
    total_impact = []
    for x in gdc_impact:
        if x[0] in cosmic_mutation:
            k = x[0]
            y = [x[0], x[1], x[2], x[3], cosmic_mutation[k][9], cosmic_mutation[k][10]]
            total_impact.append(y)
    # will only insert impact attributes for mutations found in cosmic AND gdc
    # will need to insert gdc separately
    return  total_impact

# *************************************************************************************
########## main ############

if __name__ == "__main__":

    gdc_att    = parseGDC('KIF11', 'mutations.2019-01-23.json')
    # gdc_mut    = gdc_att[0]
    # gdc_source = gdc_att[1]
    # gdc_impact = gdc_att[2]
    # #print(gdc_att[0])

    #cosmic_dict = cosmicParser('KIF11', 'V87_38_MUTANT.csv')
    #t_impact = combineImpact('KIF11', 'mutations.2019-01-23.json', 'V87_38_MUTANT.csv')
    #print(t_impact)

    i = 0
    j = 0
    for x in gdc_att[0]:
       # if x[5] != "None":  # and x[4] == "missense_variant":
            #i += 1
        print(x[4], x[3])
    #for x in cosmic_dict:
    #    j += 1
        #print(x, cosmic_dict[x])
    print('The total number of mutations in GDC is: ', i, 'and total in cosmic is: ', j)




