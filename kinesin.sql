SHOW WARNINGS;
SET SESSION sql_mode='STRICT_TRANS_TABLES';

-- CREATE DATABASE							kinesin;
USE kinesin;

DROP TABLE IF EXISTS					tissue;
DROP TABLE IF EXISTS					impact;
DROP TABLE IF EXISTS					frequency;
DROP TABLE IF EXISTS					source_info;
DROP TABLE IF EXISTS					mutation;
DROP TABLE IF EXISTS					gene;
DROP TABLE IF EXISTS					kinesin_family;

DROP VIEW IF EXISTS						vKif11_coding_mut;


CREATE table 									kinesin_family
(			family_name							VARCHAR(100)				NOT NULL,
			PRIMARY KEY (family_name)
)ENGINE=InnoDB;


CREATE table										gene
(			gene_name								VARCHAR(100)				NOT NULL,
			strand										CHAR(1)							NOT NULL,
            chrom_build							VARCHAR(10)				NOT NULL,
            chrom_location						VARCHAR(50)				NOT NULL,
            cds_length								VARCHAR(50)				NOT NULL,
			hgnc_accession						VARCHAR(50)				DEFAULT 'N/A'
																									NOT NULL,
			entrez_id									VARCHAR(50)				DEFAULT 'N/A'
																									NOT NULL,
			ensembl_id								VARCHAR(50)				DEFAULT 'N/A'
																									NOT NULL,
			synonyms								VARCHAR(200)				DEFAULT 'None'
																									NOT NULL,
			uniprot_id								VARCHAR(11)				DEFAULT 'N/A'			
																									NOT NULL,
            organism									VARCHAR(100)				NOT NULL,
            homologs								TEXT								NOT NULL,
            kinesin_family							VARCHAR(50)				NOT NULL,
            PRIMARY KEY (gene_name, organism),
            FOREIGN KEY (kinesin_family) REFERENCES kinesin_family (family_name) ON DELETE CASCADE
)ENGINE=InnoDB;


CREATE table										mutation
(			genomic_id								VARCHAR(45)				NOT NULL,
            coding										CHAR(1)							NOT NULL,
            cds											VARCHAR(10)				DEFAULT 'N/A'
																									NOT NULL,
			mutation_type							VARCHAR(20)				DEFAULT 'UNK'
																									NOT NULL,
			consequence							VARCHAR(20)				DEFAULT 'UNK'
																									NOT NULL,
            protein										VARCHAR(10)				DEFAULT 'N/A'
																									NOT NULL,
			gene_name								VARCHAR(100)				NOT NULL,
            organism									VARCHAR(100)				NOT NULL,
            domain									VARCHAR(100)				DEFAULT 'UNK'
																									NOT NULL,
            PRIMARY KEY (genomic_id),
            FOREIGN KEY (gene_name, organism) REFERENCES gene (gene_name, organism) ON DELETE CASCADE
)ENGINE=InnoDB;
           
           
CREATE table										source_info
(			source_id								VARCHAR(25)				NOT NULL,
			transcript_id							VARCHAR(25)				DEFAULT 'N/A'
																									NOT NULL,
            mutation_id								VARCHAR(45)				NOT NULL,
			PRIMARY KEY (source_id),
            FOREIGN KEY (mutation_id) REFERENCES mutation (genomic_id) ON DELETE CASCADE
)ENGINE=InnoDB;

-- deleted description table and included in mutation table
-- CREATE table										mutation_desc
-- (			mutation_id								CHAR(10)						NOT NULL,
-- 			coding 									CHAR(20)						DEFAULT 'N'
-- 																									NOT NULL,
-- 			mutation_type							VARCHAR(10)				DEFAULT 'UNK'
-- 																									NOT NULL,
-- 			consequence							VARCHAR(50)				DEFAULT 'UNK'
-- 																									NOT NULL,
-- 			PRIMARY KEY (mutation_id),
--            FOREIGN KEY (mutation_id) REFERENCES mutation (mutation_id) ON DELETE CASCADE
-- )ENGINE=InnoDB;


CREATE table										frequency
(			mutation_id								VARCHAR(45)				NOT NULL,
            thou_genome_maf  				VARCHAR(10)				DEFAULT 'UNK'
																									NOT NULL,
			goesp_maf								VARCHAR(10)				DEFAULT 'UNK'
																									NOT NULL,
			exac_maf								VARCHAR(10)				DEFAULT 'UNK'
																									NOT NULL,
			gdc_freq									VARCHAR(10)				DEFAULT 'UNK'
																									NOT NULL,
			PRIMARY KEY (mutation_id),
            FOREIGN KEY (mutation_id) REFERENCES mutation (genomic_id) ON DELETE CASCADE
)ENGINE=InnoDB;


CREATE table										impact
(			mutation_id								VARCHAR(45)				NOT NULL,
			custom_score							CHAR(10)						DEFAULT 'UNK'
																									NOT NULL,
			vep											VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			sift_score									VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			sift_prediction							VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			polyphen_score						VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			polyphen_prediction				VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			fathhm_score							VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			fathhm_prediction					VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			clinvar_prediction					VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			PRIMARY KEY (mutation_id),
            FOREIGN KEY (mutation_id) REFERENCES mutation (genomic_id) ON DELETE CASCADE
)ENGINE=InnoDB;

CREATE table										tissue
(			mutation_id								VARCHAR(45)				NOT NULL,
			tissue_type								VARCHAR(100)				DEFAULT 'UNK'
																									NOT NULL,
			cancer_type							VARCHAR(100)				DEFAULT 'UNK'
																									NOT NULL,
			PRIMARY KEY (mutation_id),
            FOREIGN KEY (mutation_id) REFERENCES mutation (genomic_id) ON DELETE CASCADE
)ENGINE=InnoDB;


-- materialised view created for user interface: includes only coding mutations for given gene (kinesin)
-- in future can expand to create new views for each gene (kinesin) or for non-coding mutations (CNVs)
CREATE VIEW vKif11_coding_mut AS 
			SELECT 	g.gene_name, g.organism, m.genomic_id, cds, protein, mutation_type, consequence, domain, custom_score
            FROM 	gene g, mutation m, impact i
            WHERE 	m.coding = 'Y'
            AND		g.gene_name = m.gene_name
			AND		g.organism = m.organism
 			AND		m.genomic_id = i.mutation_id;
            

