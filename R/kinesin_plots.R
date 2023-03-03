## Make some plots showing tissue distribution, 'consequence', domain, source info



#install.packages("ggplot2")
library(RMySQL)
library(plyr)
library(ggplot2)
library(dplyr)

### Mutations by source
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

dbListTables(mydb)
dbListFields(mydb, 'source_info')
rs <- dbSendQuery(mydb, "SELECT mutation_id source_db from kinesin.source_info
                  GROUP by source_db;")
db_distr <- dbFetch(rs, n=-1)
db_distr
#clear result set
dbClearResult(rs)

rs <- dbSendQuery(mydb, "SELECT SQL_CALC_FOUND_ROWS
                  mutation_id
                  FROM source_info 
                  GROUP by mutation_id
                  HAVING COUNT(*) > 1;")
dbClearResult(rs)
rs <- dbSendQuery(mydb, "SELECT FOUND_ROWS();")
both <- dbFetch(rs, n=-1)
dbClearResult(rs)

both
source_db <- c('BOTH')

both_db <- 52
db_sources <- c(165/257, 144/257)
install.packages('VennDiagram')
library(VennDiagram)
draw.pairwise.venn(165, 144, cross.area = 52, category = c('COSMIC', 'GDC'),
                    fill = c('blue','red'), lty='blank', cat.default.pos = 'text',
                    cat.pos = c(200, 200))

grid.newpage()
## Mutations by consequence typ
mydb=mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)
rs <- dbSendQuery(mydb, "SELECT consequence, COUNT(*) FROM kinesin.mutation 
                  GROUP by consequence;")
cons_muts <- dbFetch(rs, n=-1)
#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)

cons_muts
missense <- cons_muts[3,2]
nonsense <- cons_muts[4,2]
synonymous <- cons_muts[7,2]
splice_var <- sum(cons_muts[5:6,2])
frameshift <- cons_muts[1,2]
inframe<- cons_muts[2,2]
total <- sum(cons_muts[1:7, 2])
pie_slices <- c(missense, nonsense, synonymous, splice_var, frameshift, inframe)
pie_slices
#pie chart
cons=c('missense', 'nonsense', 'synonymous', 'splice variants', 'frameshift', 'other')
pie(pie_slices, labels=cons, col = c("purple", "violetred1", "green3", "cornsilk",
                "cyan", "white"), main='Distribution of consequence type')

    

##Mutations by domain and domain length

#fetching some data (domain breakdown)
# connect to database
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

dbListTables(mydb)
dbListFields(mydb, 'mutation')

rs <- dbSendQuery(mydb, "SELECT domain, COUNT(*) as total FROM kinesin.mutation 
                  GROUP by domain;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
domain_muts <- dbFetch(rs, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)

domain_muts
# plots of mutation numbers by domain and domain lengths
names<-c('other', 'motor', 'tail-bind')
row.names(domain_muts) <- names

domain_muts <- select(domain_muts, -1)  #remove first column
domain_muts
par_original<- par()
par(pin=c(3,3),font=2,ps=12, family="sans", cex.main=0.8, cex.lab=0.8, cex.axis=.8, mai=c(1.02, 0.82, 0.82, 0.42))
par(mfrow=c(1,2))

dom_matrix <- data.matrix(domain_muts)
dom_matrix
domain_bar<-barplot(dom_matrix, space=1, width=0.5, beside=TRUE, names.arg = names, col=rainbow(25),
               main = 'Total mutations by domain region', ylab = 'number of mutations',
               ylim = c(0,170), axes = FALSE)
axis(side = 2, at = seq(0, 200, 10), pos = 0.30)
# add x-axis with offset positions, with ticks.
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)

dom_lens <- c(582, 336, 138)
len_matrix <- data.matrix(dom_lens)

length_bar <- barplot(len_matrix, space=1.0, width=0.5, beside=TRUE, names.arg = names, col='blue', 
                      density= 10, main = 'Length of domain regions in KIF11', ylab = 'number of residues',
                      ylim = c(0, 700), axes=FALSE)
axis(side = 2, at = seq(0, 700, 50), pos = 0.30)
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)

## plot bars from mutations/length on same graph
mut_nos <- c(155, 76, 26)
mis_nos <- c(106, 56, 17)
dom_mat <- rbind(mis_nos, mut_nos)
colnames(dom_mat)<-c('other', 'motor', 'tail-bind')
dom_mat

par(mfrow=c(1,1))
mix_bar <- barplot(dom_mat, beside=TRUE, main = "Mutation distribution by domain", cex.main=1, 
                   col=c('red', 'green'), axes=FALSE, cex.names=0.6,
                   ylab = "No. of Mutations", ylim=c(0,180))
legend("topright",c("Missense","Total"), cex=0.7, 
                  fill = c("red","green"))
axis(side=2, at=seq(0,180, 20), pos=0.30, cex.axis=0.6)
axis(side=1, pos=0, tick=FALSE, labels=FALSE)

# mutation rate: ratio of total mutations/ domain length (includes synonymous mutations)
par(mfrow=c(1,1))
rat_plot <- plot(dom_lens, domain_muts[,1], type='o', xlab='Number of residues per domain', 
                 ylab='Number of mutations per domain', 
                main='Number of mutations proportional to domain length')





#binomial test to see if any enrichment in domains vs expected
# edriver function from script
eDriver <- function (MR, TM, LR, LP) {
  result <- binom.test (MR, TM, LR/LP, alternative = "greater")
  return (result$p.value)
}
sum(dom_mat[2,])
dom_mat[2,]

p_vals <- NULL
for (i in 1:length(dom_lens)) {
    mr        <- dom_mat[2,i]
    mt        <- sum(dom_mat[2,])
    p_mr_mt   <- eDriver(mr, mt, dom_lens[i], sum(dom_lens))
    p_vals    <- c(p_vals, p_mr_mt)
}
p_vals
#> p_vals
#other        motor      mt-bind 
#[1] 0.05295104 0.79888145 0.93690724
plot(dom_lens, p_vals, type='o')

## pie
par(mfrow=c(1,2))
pie(dom_matrix[,1], main='Number of mutations by domain region', labels=c('other', 'motor', 'tail-bind'))
pie(len_matrix[,1], main='Relative length of domains of KIF11', labels=c('other', 'motor', 'tail-bind'))

# graph positions of all mutations (x axis=residue position/y axis 0-5 mutations, color by tissue)
# connect to database
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

dbListTables(mydb)
dbListFields(mydb, 'mutation')

rs <- dbSendQuery(mydb, "SELECT resnum, COUNT(*), tissue_type  FROM mutation, tissue
      WHERE protein = mutation_id GROUP by resnum;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
mut_freq <- dbFetch(rs, n=-1)
mut_freq
#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)


par(mfrow=c(1,1))
freq_matrix <- data.matrix(mut_freq[,1:2])
freq_matrix
## frequency plot of mutations by residue position
freq_plot <-plot(freq_matrix, type='h',main = 'Frequency/position of mutations in KIF11', 
                ylab = 'number of mutations', ylim = c(0, 4), xlab = 'residue position', xlim = c(0, 1060))

## plot number of mutations/tissue_type
## will require categorising (have done this in eDriver_analysis.R)
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

# large intestine/colon/rectum
colres <- dbSendQuery(mydb, "SELECT resnum, COUNT(*) as colon FROM tissue t, mutation m
                      WHERE t.tissue_type in ('large_intestine', 'Colon', 'Rectum')
                      AND m.protein = t.mutation_id GROUP by resnum;")
colon_muts <- dbFetch(colres, n=-1)
colon_muts
# plot mutations/residue number
colon_matrix <- data.matrix(colon_muts)
colon_freq <- plot(colon_matrix, type='h',main = 'Frequency/position of mutations in colon', 
                  ylab = 'number of mutations', ylim = c(0, 3), xlab = 'residue position', xlim = c(0, 1060))
#add lines to indicate boundaries of motor domain
abline(v=c(24, 359), col='red')
text(180, 2.5, labels=c('motor domain'), col='red', cex=0.8)

#barplot to show distribution of missense mutations in colon tissue
rs <- dbSendQuery(mydb, "SELECT domain, COUNT(*) as total FROM mutation, tissue
                  WHERE tissue_type in ('large_intestine', 'Colon', 'Rectum')
                  AND consequence='missense'
                  AND protein=mutation_id
                  GROUP by domain;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
colon_doms <- dbFetch(rs, n=-1)
row.names(colon_doms) <- names
colon_doms <- select(colon_doms, -1)  #remove first column
colon_doms
colon_dom_matrix<-data.matrix(colon_doms)
colon_dom_matrix
# fetch results of total missense mutations by domain
rs <- dbSendQuery(mydb, "SELECT domain, COUNT(*) as total FROM kinesin.mutation 
                  WHERE consequence='missense' GROUP by domain;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
domain_missense <- dbFetch(rs, n=-1)
row.names(domain_missense) <- names
domain_missense<- select(domain_missense, -1)
dom_missense<-data.matrix(domain_missense)

#plot domain categories and number of missense mutations
par(mfrow=c(1,2))
colon_mis_bar <- barplot(colon_dom_matrix, space=1.0, width=0.5, beside=TRUE, cex.names= 0.6, 
                     names.arg = names, col='blue', density= 10, main = 'Missense Colorectal Mutations',
                     cex.main=0.8,
                     ylab = 'number of mutations', ylim = c(0, 30), axes=FALSE)
axis(side = 2, at = seq(0, 30, 5), pos = 0.30, cex.axis=0.8)
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)

domain_mis_bar <-barplot(dom_missense, space=1, width=0.5, beside=TRUE, cex.names=0.5, names.arg = names, col=rainbow(25),
                    main = 'Total missense by domain', ylab = 'missense mutations',
                    ylim = c(0,120), axes = FALSE)
axis(side = 2, at = seq(0, 200, 10), pos = 0.30)
# add x-axis with offset positions, with ticks.
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)

## plot of no missense v domain length
par(mfrow=c(1,1))
rat_plot <- plot(dom_lens, dom_missense[,1], 
                 col='blue', axes=FALSE,
                 ylim = c(0, 120), xlim=c(100, 600),
                 xlab='Number of residues per domain', 
                 ylab='Number of missense mutations per domain', 
                 main='Relationship between number of missense mutations and domain length', cex.main=0.8)
axis(side=2, at = seq(0, 100, 20), pos=90, cex.axis=0.8)
axis(side=1, pos=0, tick=TRUE, cex.axis=0.8)
text(130,30, labels="tail-binding", cex=0.6, pos=1, col="red")
text(330, 69, labels='motor', cex=0.6, pos=1, col='red')
text(590, 105, labels='other', cex=0.6, pos=1, col='red')

#correlation test
cor(dom_lens, dom_missense[,1], method='kendall')
#[1] 1
cor(dom_lens, dom_missense[,1], method='spearman')
#[1] 1

lm(dom_missense[,1] ~ dom_lens)
#Coefficients:
#(Intercept)     dom_lens  
#-10.9311       0.2006  
abline(-10.9311, 0.2006)


colon_matrix
## with all mutations total/colon
rs <- dbSendQuery(mydb, "SELECT domain, COUNT(*) as total FROM mutation, tissue
                  WHERE tissue_type in ('large_intestine', 'Colon', 'Rectum')
                  AND protein=mutation_id
                  GROUP by domain;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
colon_total <- dbFetch(rs, n=-1)
row.names(colon_total) <- names
colon_total <- select(colon_total, -1)  #remove first column
colon_total
colon_total_matrix<-data.matrix(colon_total)
colon_total_matrix
par(mfrow=c(1,2))
colon_total_bar <- barplot(colon_total_matrix, space=1.0, width=0.5, beside=TRUE, cex.names= 0.6, 
                     names.arg = names, col='blue', density= 10, main = 'All Colorectal Mutations', 
                     cex.main=0.8,
                     ylab = 'number of mutations', ylim = c(0, 30), axes=FALSE)
axis(side = 2, at = seq(0, 30, 5), pos = 0.30, cex.axis=0.8)
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)

par(mfrow=c(1,1))
domain_bar<-barplot(dom_matrix, space=1, width=0.5, beside=TRUE, cex.names=0.6, names.arg = names, col=rainbow(25),
                    main = 'Total mutations by domain', cex.main=0.8, ylab = 'mutations',
                    axes = FALSE)
axis(side = 2, at = seq(0, 200, 10), pos = 0.30, cex.axis=0.6)
# add x-axis with offset positions, with ticks.
axis(side = 1, pos = 0, tick = TRUE, labels=FALSE)


colon_dom_matrix
#plot length of domain vs number of mutations for colon missense
col_rat_plot <- plot(dom_lens, colon_dom_matrix[,1],
                 col='blue', axes=FALSE,
                 ylim = c(0, 25), xlim=c(100, 600),
                 xlab='Number of residues per domain', 
                 ylab='Number of missense mutations per domain', 
                 main='Relationship between number of missense mutations and domain length in colorectal tissue', cex.main=0.7)
axis(side=2, at = seq(0, 25, 5), pos=90, cex.axis=0.8)
axis(side=1, pos=0, tick=TRUE, cex.axis=0.8)
text(130, 5, labels="tail-binding", cex=0.6, pos=1, col="red")
text(330, 21, labels='motor', cex=0.6, pos=1, col='red')
text(590, 15, labels='other', cex=0.6, pos=1, col='red')

lm(colon_dom_matrix[,1] ~ dom_lens)
#Coefficients:
#(Intercept)     dom_lens  
#3.39559      0.02066  
abline(3.39559, 0.02066, col='black')

## correlation
cor(dom_lens, colon_dom_matrix[,1], method="spearman")
# [1] 0.5
cor(dom_lens, colon_dom_matrix[,1], method="kendall")
# [1] 0.3333333


#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)

## analyse dN/dS to see selectivity of mutations in overall gene, by domain, and by tissue and domain
# dN includes all non-synonymous mutations? not indels, so include nonsense and missense
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

rs <- dbSendQuery(mydb, "SELECT COUNT(*) FROM mutation
                  WHERE consequence = 'synonymous';")
syn<-dbFetch(rs, n=-1)
syn


ds<- dbSendQuery(mydb, "SELECT COUNT(*) FROM mutation
                 WHERE consequence = 'synonymous'
                 GROUP by domain;")
syn_dom<-dbFetch(ds, n=-1)
syn_dom
ns<- dbSendQuery(mydb, "SELECT COUNT(*) FROM mutation
                 WHERE consequence in ('missense', 'nonsense')
                 GROUP by domain;")
non_dom<-dbFetch(ns, n=-1)
non_dom
dNds1<- cbind(syn_dom, non_dom)


dNds1
dnds <- function(dN, dS) {
  result <- dN/dS
  return(round(result, 2))
}
ratio <- NULL
for (i in 1:length(dNds1)) {
  ratio <- dnds(dNds1[,2], dNds1[,1])
}
dNds <- cbind(dNds1, ratio)
colnames(dNds)<-c('synon', 'non-syn', 'ratio')

## dN/dS for entire gene: 
total_dnds <- c(sum(dNds[,1]), sum(dNds[,2]), round(sum(dNds[,2])/sum(dNds[,1]),2))
total_dnds

dNdS<- rbind(dNds, total_dnds)
dNdS
rownames(dNdS)<-c('other', 'motor', 'tail-bind', 'total')
dNdS

#         synon non-syn ratio
# other      32     112  3.50
# motor      11      61  5.55
# mt-bind     3      20  6.67
# total      46     193  4.20

# model too simplistic for comparison between genes (see Martincorena et al, 2017),
# but may be useful parameter for comparison between domains or tissues of same gene


ds<- dbSendQuery(mydb, "SELECT COUNT(*) FROM mutation, tissue
                 WHERE consequence = 'synonymous'
                 AND tissue_type in ('large_intestine', 'Colon', 'Rectum')
                 AND protein = mutation_id
                 GROUP by domain;")
syn_col<-dbFetch(ds, n=-1)
syn_col

ns<- dbSendQuery(mydb, "SELECT COUNT(*) FROM mutation, tissue
                 WHERE consequence in ('nonsense', 'missense')
                 AND tissue_type in ('large_intestine', 'Colon', 'Rectum')
                 AND protein = mutation_id
                 GROUP by domain;")
non_col<-dbFetch(ns, n=-1)
non_col


total_col <-c(sum(syn_col[,1]), sum(non_col[,1]), round(sum(non_col[,1])/sum(syn_col[,1])))
total_col
# analysis by colorectal mutations only
dNds_col <- cbind(syn_col, non_col)
dNds_col <- rbind(dNds_col, total_col)
colnames(dNds_col)<-c('synon', 'non-syn')
rownames(dNds_col)<-c('other', 'motor', 'tail-bind', 'total')
dNds_col
ratio <- NULL
for (i in 1:length(dNds_col)) {
  ratio <- dnds(dNds_col[,2], dNds_col[,1])
}
ratio
dNds_col <- cbind(dNds_col, ratio)
dNdS[,3]
ratio
ratio_matrix <- cbind(dNdS[,3], ratio)
colnames(ratio_matrix) <- c('all', 'colon')
rownames(ratio_matrix) <- c('other', 'motor', 'tail-bind', 'total')
ratio_matrix

par(mfrow=c(1,1))
ratio_plot <- plot(ratio_matrix[1:3,2], ratio_matrix[1:3,1], xlim=c(0,10), ylim=c(0,10), xlab=c('dNdS ratio colon'),
     ylab=c('dNdS ratio overall'))
lm(ratio_matrix[1:3,1] ~ ratio_matrix[1:3,2])
#Coefficients:
#(Intercept)  ratio_matrix[, 2]  
#4.2214                0.2697 
abline(4.2214, 0.2697, col='black')

text(ratio_matrix[1:3,2], ratio_matrix[1:3,1], labels = row.names(ratio_matrix[1:4,]), pos=4, cex=0.8)
title(main=c('dNdS ratio All tissues vs colon'), sub=c('(non-bootstrapped)'), cex.sub=0.6)

cor(ratio_matrix[1:3,2], ratio_matrix[1:3,1], method='kendall')
#[1] 0.3333333

#         synon non-syn ratio
# other       7      13  1.86
# motor       2      18  9.00
# mt-bind     1       2  2.00
# total      10      33  3.30

## this might be interesting because tho colorectal is known to have higher mutation rate,
#   the dN/dS ratio of mutation overall is actually lower than for the protein as a whole
#   (3.30 vs 4.20) and the difference within the regions is even more marked (ratio greater
#   within motor region).
ratio_matrix[1:3,1]

## perform kolgomorov-smirnov test to see if the ratios of colon vs overall are from same distribution
ks.test(ratio_matrix[1:3,2], ratio_matrix[1:3,1])
# D = 0.66667, p-value = 0.6
# can't reject null, may be from same distribution

## can i use wilcoxon signed ranks test in this instance (paired values)?
wilcox.test(ratio_matrix[1:3,2], ratio_matrix[1:3,1], paired=TRUE)
# V = 1, p-value = 0.5
# can't reject null, may come from same distribution

#mann whitney u test--non-paired, must be independent (NOT right test)
wilcox.test(ratio_matrix[1:3,2], ratio_matrix[1:3,1], paired=FALSE)
# W = 2, p-value = 0.4






## this is just jibberish here
ratio_matrix[,1]
ratio_matrix[,2]

## bootstrap ratios for each domain 
syn_col<-sum(dNds_col[1:3,1])
non_col<-sum(dNds_col[1:3,2])
syn_tot<-sum(dNds[,1])
syn_tot<-sum(dNds[,2])

syn_col
non_col
median_tot<-NULL
median_col<-NULL
par(mfrow=c(1,1))
for (i in 1:100){
  boot_tot<-sample(ratio_matrix[1:3,1], replace=TRUE)
  boot_col<-sample(ratio_matrix[1:3,2], replace=TRUE)
  median_tot <- c(median_tot, median(boot_tot))
  median_col <- c(median_col, median(boot_col))
  difference <- c(difference, median(boot_tot)-median(boot_col))
}
boot_col
motor_tot
median(boot_tot)
q1<-c(quantile(difference, 0.025))
q2<-c(quantile(difference, 0.975))

hist(difference, breaks=100)
abline(v=q1, col='red')
abline(v=q2, col='red')
# 0 is included in conf interval, so Ho not rejected.

plot(boot_tot, boot_col, xlim=c(0,10), ylim=c(0,10))
boot_tot

# can i bootstrap?
# test for differences between distributions in partic. tissue vs overall (total dN/ total dS)
## must be independent samples--comparing to the overall distribution not independent?

# test for differences between distributions in each domain vs overall (domain dN/ domain dS)
# test for differences between distributions in each domain of tissue vs domain of overall
# (domain-tissue DN/ domain-tissue dS)
## 1) bootstrap with replacement to create many sets of samples for each tissue vs total
## 2) compute dNdS ratio for each pair of samples (using bootstrapped samples)
## 3) compute difference of dNdS for each pair of samples
## 4) use bootstrap distribution to determine 0.95 confidence interval for dN/dS difference


## find ratios for each tissue and bootstrap overall ratios (10 values)
## not sure how this will work with too many categories.

## bootstrap values for domains?


values<-c(3.50, 1.86)
mean(values)
sd(values)

sd_ratios <- NULL
for (i in 1:4){
  sd_reg<-sd(ratio_matrix[i,])
  sd_ratios<-c(sd_ratios, sd_reg)
}
sd_ratios
sd_matrix<-matrix(sd_ratios, nrow=4, ncol=1)
sd_matrix
ratio_matrix1<-cbind(ratio_matrix, sd_ratios)
ratio_matrix1
