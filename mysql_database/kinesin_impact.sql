SHOW WARNINGS;
SET SESSION sql_mode='STRICT_TRANS_TABLES';

USE kinesin;

DROP TABLE IF EXISTS					impact;


CREATE table										impact
(			mutation_id							VARCHAR(45)				NOT NULL,
			sift_score							VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			sift_pred							VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			polyphen_score						VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			polyphen_pred       				VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
		    condel                              VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
		    cadd_rank                           VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			cadd_score                          VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
		    fathmm_rank                         VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			fathmm_score						VARCHAR(25)				DEFAULT 'UNK'
																									NOT NULL,
			fathmm_pred     					VARCHAR(25)				DEFAULT 'UNK'               NOT NULL,
			metaSVM_rank                        VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			metaSVM_score                       VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			metaSVM_pred                        VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			mutpred_rank                        VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			mutpred_score                       VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			mutassessor_rank                    VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			mutassessor_score                   VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			mutassessor_pred                    VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			muttaster_rank                      VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			muttaster_score                     VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			muttaster_pred                      VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			provean_rank                        VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			provean_score                       VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			provean_pred                        VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			revel_rank                          VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			revel_score                         VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			clinvar_prediction					VARCHAR(25)				DEFAULT 'UNK'               NOT NULL,
			median_ranked_scores                VARCHAR(25)             DEFAULT 'UNK'               NOT NULL,
			PRIMARY KEY (mutation_id),
            FOREIGN KEY (mutation_id) REFERENCES mutation (protein) ON DELETE CASCADE
)ENGINE=InnoDB;
