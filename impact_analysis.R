# Plot ranked prediction scores vs amino acid position

library(RMySQL)
library(plyr)
library(ggplot2)
library(dplyr)

#fetching some data from impact table (aa position, ranked scores)

# connect to database
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

dbListTables(mydb)
dbListFields(mydb, 'impact')

rs <- dbSendQuery(mydb, "SELECT condel, cadd_rank, fathmm_rank, metaSVM_rank
                        mutassesor_rank, muttaster_rank, provean_rank, revel_rank, resnum 
                        FROM kinesin.impact, kinesin.mutation
                        WHERE mutation_id = protein
                        ORDER  by resnum;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
impact_data <- dbFetch(rs, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)
