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

rs <- dbSendQuery(mydb, "SELECT condel, cadd_rank, fathmm_rank, metaSVM_rank,
                        mutassessor_rank, muttaster_rank, provean_rank, revel_rank, resnum 
                        FROM kinesin.impact, kinesin.mutation
                        WHERE mutation_id = protein
                        ORDER  by resnum;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
data <- dbFetch(rs, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)
View(data)
# resnum<- select(data, resnum)
# resnum
# condel<-c(select(data, condel), na.rm=TRUE)
# cadd<-c(select(data, cadd_rank))
# fathmm<-select(data, fathmm_rank)
# metaSVM<-select(data, metaSVM_rank)
# mutassessor<-select(data, mutassessor_rank)
# muttaster<-select(data, muttaster_rank)
# provean<-select(data, provean_rank)


xlim <- range(length(resnum))
y<-c(NA, 0:1, na.rm=TRUE); y
ylim <- range(y)
impact_df<-na.omit(data)
head(impact_df)

# basic scatterplot
ggplot(mpg, aes(manufacturer, cty)) + geom_point()

par(mfrow=c(1,1))
## scatter plot of various ranked scores by residue number
par_original<-par()


par_original
plot(NULL, xlim=c(0, 1060), ylim=c(0,1), ylab='ranked scores', xlab='residue',
     main='Comparison of Functional Impact Tools')
points(x=impact_df[[9]], y=impact_df[[1]], cex=0.5, col=1) #condel=black
points(x=impact_df[[9]], y=impact_df[[2]], pch=1, cex=.5, col=2) #cadd= red
#points(x=impact_df[[9]], y=impact_df[[15]], cex=0.5, col=3) #fathmm = yellow
points(x=impact_df[[9]], y=impact_df[[4]], cex=0.5, col=4) #metaSVM = blue
#points(x=impact_df[[9]], y=impact_df[[5]], cex=0.5, col=5) #mutassessor = cyan
#points(x=impact_df[[9]], y=impact_df[[6]], cex=0.5, col=6) # mutaster = pink
#points(x=impact_df[[9]], y=impact_df[[7]], cex=0.5, col=8) # provean = tan
points(x=impact_df[[9]], y=impact_df[[8]], cex=0.5, col=3) # revel = green
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
legend("topright", inset=c(-.1,0), legend=c('Condel', 'CADD', 'MetaSVM', 'Revel'), 
       pch=1, col=c(1,2,4,3), title="Tool", cex=0.6)

identify(x, label=x) # identify points 

impact_df[[1]] ## gives list of values in column 1

abline(v=impact_df[[9]])

## convert dataframe to matrix (so points are numeric)
df2 <- data.matrix(impact_df)
df2[1,9] #value of col 9, row 1

median(df2[1,1:8]) # median row 1, cols 1-8

#plot using data matrix instead of dataframe
plot(NULL, xlim=c(0, 1060), ylim=c(0,1), ylab='ranked scores', xlab='residue')
points(x=df2[,9], y=df2[,1], col=1) #condel=black

#plot for motor domain only
plot(NULL, xlim=c(0, 450), ylim=c(0,1), ylab='ranked scores', 
     xlab='residue')
points(x=df2[1:108,9], y=df2[1:108,1], col=1) #condel=black
points(x=df2[1:108,9], y=df2[1:108,2], pch=1, col=2) #cadd= red
points(x=df2[1:108,9], y=df2[1:108,3], pch=1, col=3) #fathmm = green
points(x=df2[1:108,9], y=df2[1:108,4], pch=1, col=4) #metaSVM = blue
points(x=df2[1:108,9], y=df2[1:108,5], pch=1, col=5) #mutassessor = cyan
points(x=df2[1:108,9], y=df2[1:108,6], pch=1, col=6) # mutaster = pink
points(x=df2[1:108,9], y=df2[1:108,7], pch=1, col=8) # provean = tan
points(x=df2[1:108,9], y=df2[1:108,8], pch=1, col=15) # revel = yellow
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
legend("topright", inset=c(-.15,0), legend=colnames(impact_df[-9]), pch=1, 
       col=c(1,2,3,4,5,6,8,15), title="Tool", cex=0.5)

df2[9] #gives value of col1, row 9
df2[[9]] #same

res_median<-NULL
col_median<-NULL
for (i in 1:239){
  res_median<-median(df2[i, 1:8], na.rm=TRUE)  #skips na's
  col_median<-c(col_median, res_median)
}
col_median
# plot median ranked scores 
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=FALSE)
#plot(NULL, xlim=c(0, 1060), ylim=c(0,1), ylab='ranked scores', xlab='residue', 
#     main='Ranked scores median for each residue')
plot(NULL, xaxt='n', xlim=c(0, 1060), ylim=c(0,1), xaxs='i', ylab='ranked scores', xlab='residue', 
     main='Ranked scores median for each residue')
points(x=df2[,9], y=col_median)
axis(side = 1, at = seq(0, 1060, by = 50), tick= T, labels = T, cex.axis=.6)
grid(col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
#add lines to indicate boundaries of motor domain
abline(v=c(24, 359), col='red')
text(180, 1.0, labels=c('motor domain'), col='red', cex=0.6)
abline(v=c(916, 1053), col='green')
text(985, 1.0, labels=c('mt-binding'), col='green', cex=0.6)
grid(nx = 20, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
# plot median ranked scores for motor domain only
par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=FALSE)
plot(NULL, xlim=c(0, 400), ylim=c(0,1), ylab='ranked scores', xlab='residue')
points(x=df2[1:108,9], y=col_median[1:108])

## compare sift/polyphen and condel top scores
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

rs <- dbSendQuery(mydb, "SELECT mutation_id, condel  FROM kinesin.impact
                  WHERE sift_pred='probably_damaging'
                  AND polyphen_pred = 'deleterious'
                  ORDER by condel DESC;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
sift_pp2 <- dbFetch(rs, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)
View(sift_pp2)
# 77 mutations with 'prob_damaging' and 'deleterious' predictions, 
# condel only looks at missense (56 mutations)
# what is cuttoff?
length(sift_pp2[,2])

condel_scores<-NULL
for (i in 1:77){
  condel_scores<-c(condel_scores, na.omit(sift_pp2[i,2]))
    }

min(condel_scores)  [1] "0.708"
median(condel_scores) [1] "0.945"

## analyse cadd_rank scores
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

rs <- dbSendQuery(mydb, "SELECT mutation_id, cadd_rank  FROM kinesin.impact
                  ORDER by cadd_rank DESC;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
cadd_df <- dbFetch(rs, n=-1)

#clear result set
dbClearResult(rs)
#disconnect from database
dbDisconnect(mydb)

View(cadd_df)

cadd<-NULL
cadd<-as.numeric(cadd_df[,2])
cadd<-na.omit(cadd)

## make vectors with mutations from top 20% of ranked scores (of ensemble and cadd)
## compare and use venn diagram to show which have in common
mydb = dbConnect(MySQL(), user='root', password='password', dbname='kinesin', host='localhost', port=3306)

q1 <- dbSendQuery(mydb, "SELECT mutation_id FROM impact
                        WHERE condel > 0.8 GROUP by mutation_id;")
# fetch results ('n=-1' retrieves all pending records), saved as dataframe obj
condel_rs <- dbFetch(q1, n=-1)
dbClearResult(q1)
q2<- dbSendQuery(mydb, "SELECT mutation_id FROM impact
                        WHERE cadd_rank > 0.8 GROUP by mutation_id;")
cadd_rs<- dbFetch(q2, n=-1)
dbClearResult(q2)
q3<- dbSendQuery(mydb, "SELECT mutation_id FROM impact
                        WHERE metaSVM_rank > 0.8 GROUP by mutation_id;")
meta_rs<- dbFetch(q3, n=-1)
dbClearResult(q3)
q4<-dbSendQuery(mydb, "SELECT mutation_id FROM impact
                        WHERE revel_rank > 0.8 GROUP by mutation_id;")
revel_rs<-dbFetch(q4, n=-1)
dbClearResult(q4)
q5<-dbSendQuery(mydb, "SELECT mutation_id FROM impact
                        WHERE fathmm_cancer_pred = 'CANCER' GROUP by mutation_id;")
fathmm_rs<-dbFetch(q5, n=-1)
dbClearResult(q5)
#disconnect from database
dbDisconnect(mydb)

condel_top<-condel_rs$mutation_id
cadd_top<-cadd_rs$mutation_id
meta_top<-meta_rs$mutation_id
revel_top<-revel_rs$mutation_id
fathmm_cancer<-fathmm_rs$mutation_id

