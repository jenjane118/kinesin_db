# kinesin_db

Jennifer J. Stiens

Kinesin mutation database project

MSc in Bioinformatics, Birkbeck College, University of London

j.j.stiens@gmail.com


This is a thesis project for the MSc in Bioinformatics at Birkbeck College, University of London, supervised by Carolyn
Moores and Maya Topf. It involves parsing COSMIC and TCGA databases for somatic mutations in cancer found in Kinesin-5
of genes, namely the human gene known as 'KIF11' or 'Eg5'. The relational database is a MySQL database and the other
scripts are written in python. Statistical analysis is in R.

To create database, download and run <kinesin.sql>

Download the following files: 

<mutation_parser.py>, <dbinsert_module.py>, <impact_table.py>, <domain_mapper.py>, <median_score.py>, config file (like <config_kinesin.py>), and data files (<mutations.2019-01-23.json>, <results.json>, <V87_38_MUTANT.csv>, <clinvar_result.txt>, <fathmm_results.txt>, <vep_complete_results.txt>)

To populate database, run:

<mutation_parser.py>, <dbinsert_module.py>, <impact_table.py>, <vep_parse.py>


