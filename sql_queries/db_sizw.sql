

SELECT schema as "database name"
sum(data_length+index_length)/1024/1024 as "Database Size in MB"
from information_schema.tables
where table_schema = 'kinesin';
