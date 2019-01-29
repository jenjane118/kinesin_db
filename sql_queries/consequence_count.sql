
SELECT consequence, COUNT(*) AS count
FROM kinesin.mutation
GROUP BY consequence;
