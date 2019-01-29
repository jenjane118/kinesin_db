USE kinesin;
SELECT source_db, COUNT(*) as count
FROM source_info
GROUP BY source_db;

 

