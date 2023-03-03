USE kinesin;

-- changed preferences out of 'safe updates' mode to allow deleting all contents; may want to put back (prefs, sql editor, reconnect to db) when adding real data.

DELETE FROM kinesin_family;
DELETE FROM gene;
DELETE FROM mutation;
DELETE FROM source_info;
DELETE FROM impact;
DELETE FROM tissue;

INSERT INTO	 kinesin_family (family_name)
					VALUES ('Kinesin-5');

INSERT INTO 	gene (gene_name, strand, chrom_build, chrom_location, cds_length, hgnc_accession, entrez_id, ensembl_ID, synonyms, uniprot_id, organism, homologs, kinesin_family)
					VALUES ('KIF11', '+', 'h38', '10:92,593,286 - 92,655,395', '3171', '6388', '3832', 'ENSG00000138160', 'Eg5, HKSP, KNSL1, TRIP5', 'P52732', 'homo sapiens', 'mus musculus, rattus norwegicus', 'Kinesin-5');

INSERT INTO mutation (genomic, coding, cds, mutation_type, consequence, protein, gene_name, organism)
					VALUES ('10:92616824C>T', 'Y', '1120', 'substitution', 'missense', 'L374F', 'KIF11', 'homo sapiens');
                    
INSERT INTO source_info (source_id, source_db, transcript_id, mutation_id)
					VALUES ('COSM372619', 'COSMIC', 'ENST00000260731', 'L374F');
                    
INSERT INTO	impact (mutation_id, fathhm_score, fathhm_prediction)
					VALUES ('L374F', '0.098623', 'pathogenic');
                    
INSERT INTO tissue (mutation_id, tissue_type, cancer_type)
					VALUES ('L374F', 'lung', 'adenocarcinoma');
                   