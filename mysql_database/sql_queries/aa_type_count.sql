USE kinesin;
SELECT SUBSTR(protein, 1,1) AS alpha,  COUNT(*) AS count
FROM mutation
GROUP BY SUBSTR(protein, 1,1);


