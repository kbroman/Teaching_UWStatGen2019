## Basics of QTL mapping

## 1: load data + summaries

# load R/qtl
library(qtl)

# load example data set from web
sug <- read.cross("csv", "http://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

# summary information
summary(sug)

# other summary info
nind(sug)
nmar(sug)
totmar(sug)
nphe(sug)
nchr(sug)

# summary plot
plot(sug)

# individual parts of that plot
plotMap(sug)
plotMissing(sug)
plotPheno(sug, 1)
plotPheno(sug, 2)

## 2: calculate genotype probabilities; QTL genome scan

# calculate QTL genotype probabilities
sug <- calc.genoprob(sug, step=1)

# genome scan of first phenotype by interval mapping
out.em <- scanone(sug)

# plot LOD curves
plot(out.em)

# interactive plot
library(qtlcharts)
iplotScanone(out.em, sug, chr=c(7, 11, 15))

# top LOD score on each chromosome
summary(out.em)

# LOD scores above 4
summary(out.em, threshold=4)

## 3. permutation test
# we'll just do 200 permutations, because it takes a while
operm.em <- scanone(sug, n.perm=200)

# significance thresholds
summary(operm.em)

# alpha = 0.05, 0.2
summary(operm.em, alpha=c(0.05, 0.2))

# histogram of results
plot(operm.em)

# Add threshold to plot
plot(out.em)
add.threshold(out.em, perms=operm.em, lty=2, col="orchid")

## 4. LOD support intervals
# default is to drop 1.5
lodint(out.em, chr=7)

# 2-LOD support interval
lodint(out.em, chr=7, drop=2)

# expand to flanking markers
lodint(out.em, chr=7, drop=2, expandtomarkers=TRUE)

# chr 15
lodint(out.em, chr=15, drop=2, expandtomarkers=TRUE)

# approximate Bayes intervals
bayesint(out.em, chr=7, expandtomarkers=TRUE)
bayesint(out.em, chr=15, expandtomarkers=TRUE)

## 5. Haley-Knott regression
out.hk <- scanone(sug, method="hk")

# plot the two together
plot(out.em, out.hk, col=c("slateblue", "orchid"), lty=c(1,2))

# another way to make that plot
plot(out.em, col="slateblue")
plot(out.hk, col="orchid", lty=2, add=TRUE)

# plot the differences
plot(out.em - out.hk, ylim=c(-0.5, 0.5), ylab="LOD(EM)-LOD(HK)")

# real advantage is with permutations
operm.hk <- scanone(sug, method="hk", n.perm=1000)

# can also use multiple CPU
operm.hk <- scanone(sug, method="hk", n.perm=1000, n.cluster=8)
