SELECT t.mutation_id, m.consequence, m.mutation_type
FROM mutation m, tissue t
WHERE t.tissue_type = 'large_intestine';