
SELECT mutation_id, COUNT(*)
FROM  source_info
GROUP BY mutation_id
		HAVING COUNT(*) >1;
