SELECT mutation_id, condel, cadd_rank, metaSVM_rank, revel_rank, fathmm_cancer_pred
                  FROM kinesin.impact, kinesin.mutation
                  WHERE condel > 0.8 OR cadd_rank>0.8 OR metaSVM_rank>0.8 OR revel_rank>0.8 OR fathmm_cancer_pred = 'CANCER'
                  GROUP by mutation_id;
                  
                  
SELECT mutation_id FROM impact
WHERE condel > 0.8
GROUP by mutation_id;

SELECT mutation_id FROM impact
WHERE cadd_rank > 0.8
GROUP BY mutation_id;

SELECT mutation_id FROM impact
WHERE metaSVM_rank > 0.8
GROUP BY mutation_id;

SELECT mutation_id FROM impact
WHERE revel_rank > 0.8
GROUP BY mutation_id;

SELECT mutation_id FROM impact
WHERE fathmm_cancer_pred = 'CANCER'
GROUP BY mutation_id;

SELECT mutation_id FROM impact
WHERE fathmm_cancer_pred = 'CANCER'
AND clinvar_prediction in ('Pathogenic', 'Likely pathogenic');

SELECT mutation_id, condel, cadd_rank, fathmm_cancer_pred
                  FROM kinesin.impact, kinesin.mutation
                  WHERE  cadd_rank >0.8   AND fathmm_cancer_pred = 'CANCER'
                  GROUP by mutation_id;
                  
SELECT * from tissue
WHERE mutation_id in ('L991V', 'F1012I', 'L500F', 'P1010Q');

SELECT mutation_id, domain
                  FROM kinesin.impact, kinesin.mutation
                  WHERE condel > 0.8 AND cadd_rank>0.8 AND metaSVM_rank>0.8 AND revel_rank>0.8
                  AND mutation_id = protein
                  GROUP by mutation_id
                  ORDER by resnum ASC;
                  
SELECT * FROM impact
WHERE median_ranked_scores >=0.8;

SELECT mutation_id, domain
FROM kinesin.impact, kinesin.mutation
                  WHERE condel > 0.8 AND cadd_rank>=0.8 AND metaSVM_rank>=0.8 AND revel_rank>=0.8 AND median_ranked_scores>=0.8
                  AND mutation_id = protein
                  GROUP by mutation_id
                  ORDER by resnum ASC;
