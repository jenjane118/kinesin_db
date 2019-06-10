#!/usr/bin python3

# function to parse fathmm results file

import csv

def fathmmResultsParser(csv_file):

    fathmm_dict = {}

    with open (csv_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        row_one = True
        for row in csv_reader:
            if row_one == False:
                try:
                    mutation    = row[3]
                    prediction  = row[4]
                    score       = row[5]
                    fathmm_dict[mutation] = prediction, score
                except NameError as e:
                    print("Error", e)
            else:
                row_one = False
    file.close()

    return fathmm_dict


########## main ############

if __name__ == "__main__":

    results = fathmmResultsParser('fathmm_results.txt')
    print(results)

