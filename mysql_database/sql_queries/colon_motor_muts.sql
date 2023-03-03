## query for colorectal missense mutations in motor domain
SELECT	protein, median_ranked_scores FROM mutation m, tissue t, impact i
WHERE consequence = 'missense'
AND domain = 'kinesin motor'
AND tissue_type in ('Colon', 'large_intestine', 'Rectum')
AND m.protein = i.mutation_id
AND i.mutation_id = t.mutation_id
GROUP by resnum;