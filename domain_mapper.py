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
V1.0    29.01.10        Initial version                             By: JJS

"""

# ******************************************************************************
# Import libraries

import re
import sys

# ******************************************************************************

def domainMapper(mutation):
    """This function maps the amino acid mutation to a domain on the KIF11-Human pfam architecture.
    It uses uniprot id P52732
    Input               mutation            mutated amino acid residue
    Output              domain_name         name of relevant domain
                                            'motor', 'microtubule-binding', 'UNK'
    """
    import re

    domain_dict = {}
    for x in range(1,2000):
        if x in range(24,359):
            domain_dict[str(x)] = 'kinesin motor'
        elif x in range(916, 1053): #> 915 and x < 1054:
            domain_dict[str(x)] = 'microtubule-binding'
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

########## main ############

if __name__ == "__main__":

    p = domainMapper('A366T')
    print(p)