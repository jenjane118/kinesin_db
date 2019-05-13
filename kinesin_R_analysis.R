# Kinesin mutation analysis
# j.j.stiens@gmail.com

# 08.05.19

# Use right-sided binomial test to see if mutations are enriched within any domains
# in any particular tissues in kinesin-5.

# Run test for each tissue sample in database. Use tissue rather than cancer type
# e-driver analyses using cancer types from TCGA which are specific to tissue (ex: breast
#     sarcoma, uterine carcinoma, etc)

# use RMySQL package to query straight from R
#install.packages("RMySQL")
#library(RMySQL)
library(plyr)
install.packages('dplyr')
library(dplyr)
# connect to database
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

dbListTables(mydb)
dbListFields(mydb, 'mutation')

# make dataframe of number of mutations in each tissue, grouped by domain.
# make dataframe of tissues/number of mutations/domain
#  columns: tissue type
#  rows: domain names

# Query for total number of missense mutations in each domain
rs <- dbSendQuery(mydb, "SELECT domain, COUNT(*) as total FROM kinesin.mutation 
                  WHERE consequence = 'missense' GROUP by domain;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
total_mis <- dbFetch(rs, n=-1)
# large intestine/colon/rectum
colres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as colon FROM tissue t, mutation m
                      WHERE consequence = 'missense'
                      AND tissue_type in ('large_intestine', 'Colon', 'Rectum')
                      AND m.protein = t.mutation_id
                      GROUP BY m.domain;")
colon_mis <- dbFetch(colres, n=-1)
# Bladder/urinary tract
ut <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as urinary FROM tissue t, mutation m
                  WHERE consequence = 'missense'
                  AND tissue_type in ('urinary_tract', 'Bladder')
                  AND m.protein = t.mutation_id
                  GROUP BY m.domain;")
urinary_mis <- dbFetch(ut, n=-1)
# Corpus uteri/uterus/endometrium
cu <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as uterus FROM tissue t, mutation m
                  WHERE consequence = 'missense'
                  AND tissue_type in ('Corpus uteri', 'Uterus, NOS', 'endometrium')
                  AND m.protein = t.mutation_id
                  GROUP BY m.domain;")
uterus_mis <- dbFetch(cu, n=-1)
# Liver/biliary_tract/liver and bile ducts
livres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as liver FROM tissue t, mutation m
                      WHERE consequence = 'missense'
                      AND tissue_type in ('liver', '%liver', 'biliary_tract')
                      AND m.protein = t.mutation_id
                      GROUP BY m.domain;")
liver_mis <- dbFetch(livres, n=-1)
# upper aerodigestive tract/Tonsil/oesophagus/Esophagus
uat <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as upresp FROM tissue t, mutation m
                   WHERE consequence = 'missense'
                   AND tissue_type in ('%upper', 'Tonsil', 'oesophagus', 'Esophagus', 'Floor of mouth')
                   AND m.protein = t.mutation_id
                   GROUP BY m.domain;")
upres_mis <- dbFetch(uat, n=-1)
# prostate/prostate gland
pres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as prostate FROM tissue t, mutation m
                    WHERE consequence = 'missense'
                    AND tissue_type in ('prostate', '%Prostate')
                    AND m.protein = t.mutation_id
                    GROUP BY m.domain;")
prostate_mis <- dbFetch(pres, n=-1)
#breast
bres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as breast FROM tissue t, mutation m
                    WHERE consequence = 'missense'
                    AND tissue_type = 'breast'
                    AND m.protein = t.mutation_id
                    GROUP BY m.domain;")
breast_mis <- dbFetch(bres, n=-1)
#skin
skres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as skin FROM tissue t, mutation m
                     WHERE consequence = 'missense'
                     AND tissue_type LIKE 'skin'
                     AND m.protein = t.mutation_id
                     GROUP BY m.domain;")
skin_mis <- dbFetch(skres, n=-1)
#lung
lres <- dbSendQuery(mydb, "SELECT m.domain, COUNT(*) as lung FROM tissue t, mutation m
                    WHERE consequence = 'missense'
                    AND tissue_type LIKE 'lung'
                    AND m.protein = t.mutation_id
                    GROUP BY m.domain;")
lung_mis <- dbFetch(lres, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)

#merge data into one dataframe
tissue_mis <-NULL
tissue_mis <- merge(total_mis, colon_mis)
tissue_mis <- merge(tissue_mis, uterus_mis)
tissue_mis <-merge(tissue_mis, prostate_mis)
tissue_mis <-merge(tissue_mis, skin_mis)
tissue_mis <- merge(tissue_mis, breast_mis)
tissue_mis <- merge(tissue_mis, lung_mis)
tissue_mis <- merge(tissue_mis, urinary_mis, all.x=TRUE)    #will sub in na for missing values
tissue_mis <- merge(tissue_mis, upres_mis, all.x=TRUE)
tissue_mis <- merge(tissue_mis, liver_mis, all.x=TRUE)

row.names(tissue_mis) <- c('other', 'motor', 'mt-bind')
tissue_missense <- select(tissue_mis, -1)  #remove first column
View(tissue_missense)

######### binomial test
## binomial test to identify 'driver PFRs'
## evaluate distr of mutations within regions compared to random distr
## outputs dataframe with P-values of mutation distribution using adapted Fisher's test
## Porta-Pardo and Godzik, 2014

# edriver function from script
eDriver <- function (MR, TM, LR, LP) {
  result <- binom.test (MR, TM, LR/LP, alternative = "greater")
  return (result$p.value)
}
# run test on each tissue for each domain
dom_lens <- c(582, 336, 138)
kif11_len <- 1056
res_missense <- NULL
for (j in 1:10) {
  p_vals <- NULL
  for (i in 1:3) {
    if (is.na(tissue_missense[i,j])) {
      p_mr_mt <- NA
    }else {
      mr        <- tissue_missense[i,j]
      mt        <- sum(tissue_missense[j], na.rm=TRUE)
      p_mr_mt   <- eDriver(mr, mt, dom_lens[i], kif11_len)
    }
    p_vals    <- c(p_vals, p_mr_mt)
  }
  if (j==1) {
    res_missense <- data.frame(p_vals, row.names=c('other', 'motor', 'mt-bind'))
  }else {
    res_missense[j] <- c(p_vals)
  }
}
colnames(res_missense) <- c('total', 'colon', 'uterus', 'prostate', 'skin', 'breast', 'lung', 'upresp', 'urinary', 'liver')
View(res_missense)

## use idr1 751-885 to see if enriched (24 mutations)
idr1 <- eDriver(24, 179, 135, 1056)
idr1
# [1] 0.9999995      NOT enriched in this region.
## idr2 422-496
idr2 <- eDriver(14, 179, 74, 1056)
idr2
#[1] 0.3750199  
# use only with colon mutations (five mutations from colon)
idr1_colon <- eDriver(5, 32, 135, 1056)
idr1_colon
#[1] 0.389876   Not significant. 
## no mutations in idr2 region in colon




## test with calculations from paper

e_res <- eDriver(29, 43, 227, 732)
e_res
#[1] 9.782058e-07  doesn't match value in Figure 1
binom.test(29, 43, 227/732, alternative='greater')
#  p-value = 9.782e-07
## testing calculations in paper
mt <- 43
mr <- 29
P_mutreg <- 0.31
p_test   <- choose(mt,mr)*(P_mutreg^mr)*((1-P_mutreg)^(mt-mr))
p_test
#[1] 7.719264e-07  doesn't match value in Figure 1 (is this probability vs p_value?)



#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)
