#!/usr/bin python3

""" Impact Table Parsing and Population """

"""
Program:    impact_table
File:       impact_table.py
Version:    1.0
Date:       18.06.19
Function:   Parse VEP, FATHHM and Clinvar files and insert results into impact table of kinesin database
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program parses VEP (Variant Effect Predictor) files (https://www.ensembl.org/Homo_sapiens/Tools/VEP), FATHMM 
(http://fathmm.biocompute.org.uk) and Clinvar (https://www.ncbi.nlm.nih.gov/clinvar)for impact results. 
It includes functions to parse the output files and update mysql database impact table.

Usage:
======
impact_table            SELF

Revision History:
=================
V1.0    09.04.19        Initial: combined three scripts (vep_parse.py, fathmm_impact.py and     By: JJS
                        clinvar_impact.py) to make one module for impact table.        

"""

#*****************************************************************************
# Import libraries

import csv
import re
import pymysql
import config_home
import config_kinesin
import Bio.Data.IUPACData

#*****************************************************************************
def parseVep2(my_gene, vep_file):
    """ Parse VEP flat file results (downloaded after using webservice (https://www.ensembl.org/Multi/Tools/VEP?db=core)
    with file from vep_parse.vepScores function) to get all predictions and scores for
    entry to kinesin database impact table.

    Input               my_gene                 Gene of interest (KIF11)
                        vep_file                flat file (results) downloaded from VEP website
                                                ('vep_complete_results.txt')
    Output              impact_list             list of attributes for impact table update
    """

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
                genomic         = row['#Uploaded_variation']
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

def fathmmResultsParser(csv_file):
    """This function parses FATHMM results file downloaded from FATHMM website (http://fathmm.biocompute.org.uk)
    for predictions of the oncogenic status of missense amino-acid substitutions in the KIF11 protein. This file includes
    all missense mutations for the gene 'KIF11' and assesses 'driver' or 'passenger' status. The prediction threshold for
    the algorithm was set for '1.0', which is slightly more sensitive and less-specific than the default values.

        Input               csv_file            file of FATHMM results returned from web query
                                                ('fathmm_results.txt')
        Output              fathmm_dict         dictionary mutation identifier:prediction, score
        """

    fathmm_dict  = {}
    with open (csv_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        row_one = True
        for row in csv_reader:
            if row_one == False:
                try:
                    mutation     = row[3]
                    prediction   = row[4]
                    score        = row[5]
                    fathmm_dict[mutation] = prediction, score
                except NameError as e:
                    print("Error", e)
            else:
                row_one = False
    file.close()

    return fathmm_dict

#*****************************************************************************

def fathmmInsert(results_dict, database):
    """ This function uses list of clinvar attributes to update impact table in kinesin database.
        Input               results_dict                    dict of attributes (mutation_id:prediction, score)
                            database                        which database used (home or kenobi)
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

    sql_fathmm2 = "UPDATE impact SET fathmm_cancer_pred = %s, fathmm_cancer_score = %s WHERE mutation_id = %s;"

    count = 0
    for k,v in results_dict.items():
        aa      = k
        prediction = v[0]
        score   = v[1]
        entry   = (prediction, score, aa)
        with cnx.cursor() as cursor:
            cursor.execute(sql_fathmm2, entry)
            count += 1

    cnx.commit()
    cnx.close()

    return count

# ******************************************************************************

def parseClinvar(gene, csv_file):
    """This function parses Clinvar file (https://www.ncbi.nlm.nih.gov/clinvar) for clinical significance. This file
    includes all recorded mutations for the gene 'KIF11' that involve only a single gene. Most of these are germline
    variants and won't be present in kinesin db of somatic mutations from cancer. It includes functions to parse the
    output file from clinvar and update mysql database tables with clinical significance.

    Input               csv_file            file of Clinvar results returned from web query ('clinvar_result.txt')
    Output              cv_list             list of mutation identifier and clinical significance
    """

    cv_list  = []

    with open(csv_file, 'r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')

        for row in csv_reader:

            clin_sig = row['Clinical significance (Last reviewed)']
            # remove review date in parentheses
            r = re.compile(r'\(.*\)$')
            clin_sig = r.sub('', clin_sig)

            mut_name    = row['Name']       ##NM_004523.3(KIF11):c.2T>C (p.Met1Thr)
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

    gene = 'KIF11'
    db = 'home'

    impact = parseVep2(gene, 'vep_complete_results.txt')
    number = updateImpact2(impact, db)
    print(number)
    results = fathmmResultsParser('fathmm_results.txt')
    rows = fathmmInsert(results, db)
    print(rows)
    results2 = parseClinvar(mygene, 'clinvar_result.txt')
    successful = clinvarUpdate(results2, db)
    print(successful)