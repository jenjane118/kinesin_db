#!/usr/bin python3

""" MEDIAN ranked scores """

"""
Program:    median_score
File:       median_score.py
Version:    1.0
Date:       14.06.19
Function:   Calculate median of ranked impact scores and insert into impact table in kinesin db
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/kinesin_db

Course:     MSc Bioinformatics, Birkbeck University of London
            Project supervisors:  Dr. Carolyn Moores
                                  Dr. Maya Topf
_____________________________________________________________________________
Description:
============
This program calculates the median of the ranked scores from all the tools in the kinesin db.
 There are currently 8 ranked scores for each mutation, plus Condel, which is an ensemble score
 based on Sift and Polyphen2 which is normalised to 0-1 scale. I will use this in the median values.
 

Usage:
======
median_score                SELF

Revision History:
=================
V1.0    14.06.19        Initial                                     By: JJS

"""

#
#*****************************************************************************
# Import libraries

import re
import config_home
import pymysql
import config_kinesin
import math


#*****************************************************************************


def getScores(database):
    """Query kinesin database and get relevant ranked scores for each mutation in
    the impact table.
   INPUT            database                        home or kinesin
   OUTPUT           score_dict                      dictionary of mutation:scores
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

    mutation_list = []

    query = "SELECT mutation_id, condel, cadd_rank, fathmm_rank, metaSVM_rank, mutpred_rank, mutassessor_rank,"\
            "muttaster_rank, provean_rank, revel_rank FROM kinesin.impact;"

    with cnx.cursor() as cursor:
        cursor.execute(query)
        temp = cursor.fetchall()

    # create dictionary with scores
    mut_scores = ''
    score_dict = {}

    for x in temp:
        mut_scores = x[1:10]
        score_dict[x[0]]= mut_scores

    return score_dict


#*****************************************************************************

def calcMedian(score_dict):
    """Calculate median score for each mutation in score_dict. Median value of 9 scores.
    INPUT           score_dict                      dictionary of mutations:scores
    OUTPUT          median_dict                     dictionary of mutation:median
    """
    median_dict  = {}

    for k,v in score_dict.items():
            scores          = v[0:9]
            score_list      = []
            median_score    = ''
            for x in scores:
                if x != 'NA':
                    score_list.append(x)
                    score_list.sort()
                    median_score = score_list[math.floor(len(score_list) / 2)]
                else:
                    median_score = 'NA'
            median_dict[k] = median_score


    return median_dict



#*****************************************************************************

def insertMedians(median_dict, database):
    """Insert median values into impact database (column name: median_ranked_scores)
    INPUT           median_dict                     dictionary of mutation:median
    OUTPUT          number                          number of rows inserted
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

    sql_median = "UPDATE impact SET median_ranked_scores = %s WHERE mutation_id = %s;"

    median_list = []
    count = 0
    for k,v in median_dict.items():
        median_list   = (v,k)
        with cnx.cursor() as cursor:
            rows = cursor.execute(sql_median, median_list)
            count += 1

    cnx.commit()
    cnx.close()

    return count


# ******************************************************************************


########## main ############

if __name__ == "__main__":


    dict = getScores('home')

    medians = calcMedian(dict)
    

    rows = insertMedians(medians, 'home')
    print(rows)
