rm(list=ls())
library(extrafont)
library(showtext)
font_add_google("Roboto Slab", "roboto_slab")
showtext_auto()

pluri.axd9 <- read.csv('axd9.pluri.csv')
str(pluri.axd9)
samples <- c(rep('axd9',3), rep('h1',3))
mat.data <- cbind(samples, pluri.axd9[,2], pluri.axd9[,3], pluri.axd9[,4])
mat.data[1:6,2:4] <- as.numeric(mat.data[1:6,2:4])

col.vec <- c('samples', 'SOX2', 'NANOG', 'OCT4')
colnames(mat.data) <- col.vec
df.data <- data.frame(mat.data)
str(df.data)
df.data$samples <- factor(df.data$samples)
# df.data$GadPH <- as.numeric(df.data$GadPH)
df.data$SOX2 <- as.numeric(df.data$SOX2)
df.data$NANOG <- as.numeric(df.data$NANOG)
df.data$OCT4 <- as.numeric(df.data$OCT4)
str(df.data)
# data.G <- aggregate(list('mean.GadPH'=df.data$GadPH), by=list('Sample'=df.data$samples), FUN=mean)
data.S <- aggregate(list('mean.SOX2'=df.data$SOX2), by=list('Sample'=df.data$samples), FUN=mean)
data.N <- aggregate(list('mean.NANOG'=df.data$NANOG), by=list('Sample'=df.data$samples), FUN=mean)
data.O <- aggregate(list('mean.OCT4'=df.data$OCT4), by=list('Sample'=df.data$samples), FUN=mean)

std.err <- function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
# se.G <- aggregate(list(se.GadPH=df.data$GadPH), by=list(samples=df.data$samples), std.err)
se.S <- aggregate(list(se.SOX2=df.data$SOX2), by=list(samples=df.data$samples), std.err)
se.N <- aggregate(list(se.NANOG=df.data$NANOG), by=list(samples=df.data$samples), std.err)
se.O <- aggregate(list(se.OCT4=df.data$OCT4), by=list(samples=df.data$samples), std.err)


# GadPH <- data.G$mean.GadPH
SOX2 <- data.S$mean.SOX2
NANOG <- data.N$mean.NANOG
POU5F1 <- data.O$mean.OCT4
mat.mean <- cbind(SOX2, NANOG, POU5F1)

# se.GadPH <- se.G$se.GadPH
se.SOX2 <- se.S$se.SOX2
se.NANOG <- se.N$se.NANOG
se.OCT4 <- se.O$se.OCT4
mat.se <- cbind(se.SOX2, se.NANOG, se.OCT4)

par(family="roboto_slab")

bp<-barplot(mat.mean, beside=TRUE, ylim=c(0,0.12), xlab='Genes analysed', ylab='Expression normalized to GadPH', main='Pluripotency analysis of LU AxD9', col=c('lightseagreen', 'lightcoral'), legend.text =c('LU AxD9', 'H1 positive control'))
arrows(bp, mat.mean+mat.se, bp, mat.mean-mat.se, code=3, angle=90)
