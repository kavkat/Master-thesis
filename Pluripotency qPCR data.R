rm(list=ls())
library(extrafont)
library(showtext)
font_add_google("Roboto Slab", "roboto_slab")
showtext_auto()

#reading data
pluri.axd <- read.csv('axd.pluri.csv')
str(pluri.axd) #check data types
samples <- c(rep('axd1',3), rep('h1.1',3),rep('axd9',3), rep('h1.9',3)) 
mat.data <- cbind(samples, pluri.axd[,2], pluri.axd[,3], pluri.axd[,4])
mat.data[1:12,2:4] <- as.numeric(mat.data[1:12,2:4])

#building data frame
col.vec <- c('samples', 'SOX2', 'NANOG', 'OCT4')
colnames(mat.data) <- col.vec
df.data <- data.frame(mat.data)
str(df.data)
df.data$samples <- factor(df.data$samples)
df.data$SOX2 <- as.numeric(df.data$SOX2)
df.data$NANOG <- as.numeric(df.data$NANOG)
df.data$OCT4 <- as.numeric(df.data$OCT4)
str(df.data)

#aggregate data based on sample and calculate means
data.S <- aggregate(list('mean.SOX2'=df.data$SOX2), by=list('Sample'=df.data$samples), FUN=mean)
data.N <- aggregate(list('mean.NANOG'=df.data$NANOG), by=list('Sample'=df.data$samples), FUN=mean)
data.O <- aggregate(list('mean.OCT4'=df.data$OCT4), by=list('Sample'=df.data$samples), FUN=mean)

std.err <- function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))

#calculate standard error values for each sample
se.S <- aggregate(list(se.SOX2=df.data$SOX2), by=list(samples=df.data$samples), std.err)
se.N <- aggregate(list(se.NANOG=df.data$NANOG), by=list(samples=df.data$samples), std.err)
se.O <- aggregate(list(se.OCT4=df.data$OCT4), by=list(samples=df.data$samples), std.err)

#matrix data of means
SOX2 <- data.S$mean.SOX2
NANOG <- data.N$mean.NANOG
POU5F1 <- data.O$mean.OCT4
mat.mean <- cbind(SOX2, NANOG, POU5F1)

#matrix data of SE
se.SOX2 <- se.S$se.SOX2
se.NANOG <- se.N$se.NANOG
se.OCT4 <- se.O$se.OCT4
mat.se <- cbind(se.SOX2, se.NANOG, se.OCT4)

par(mfrow=c(1,1), family="roboto_slab")

bp<-barplot(mat.mean[,1], beside=TRUE, ylim=c(0,0.12), xlab='Samples analysed', ylab='Expression normalized to GadPH', main='Sox2 Expression data', col=c('lightseagreen','lightcoral','yellow', 'wheat'), names.arg =c('LU AxD1', 'H1 positive \ncontrol(AxD1)', 'LU AxD9', 'H1 postitive \ncontrol(AxD9)'))
arrows(bp, mat.mean[,1]+mat.se[,1], bp, mat.mean[,1]-mat.se[,1], code=3, angle=90)

bp<-barplot(mat.mean[,2], beside=TRUE, ylim=c(0,0.12), xlab='Samples analysed', ylab='Expression normalized to GadPH', main='NANOG Expression data', col=c('lightseagreen','lightcoral','yellow', 'wheat'), names.arg =c('LU AxD1', 'H1 positive \ncontrol(AxD1)', 'LU AxD9', 'H1 postitive \ncontrol(AxD9)'))
arrows(bp, mat.mean[,2]+mat.se[,2], bp, mat.mean[,2]-mat.se[,2], code=3, angle=90)

bp<-barplot(mat.mean[,3], beside=TRUE, ylim=c(0,0.009), xlab='Samples analysed', ylab='Expression normalized to GadPH', main='POU5F1 (Oct4) Expression data', col=c('lightseagreen','lightcoral','yellow', 'wheat'), names.arg =c('LU AxD1', 'H1 positive \ncontrol(AxD1)', 'LU AxD9', 'H1 postitive \ncontrol(AxD9)'))
arrows(bp, mat.mean[,3]+mat.se[,3], bp, mat.mean[,3]-mat.se[,3], code=3, angle=90)

#statistical significance - one-sample t-test
#testing normality for all data sets
##SOX2
shapiro.test(df.data$SOX2[1:3]) #not normal
shapiro.test(df.data$SOX2[7:9]) #normal
##NANOG
shapiro.test(df.data$NANOG[1:3]) #normal
shapiro.test(df.data$NANOG[7:9]) #normal
##Oct4
shapiro.test(df.data$OCT4[1:3]) #normal
shapiro.test(df.data$OCT4[7:9]) #not normal

#assumption of normality of t-test not fulfilled for all data

#use non-parametric test - Wilcoxon signed ranks test
##SOX2
wilcox.test(df.data$SOX2[1:3], mu=mat.mean[2,1])
#null-hypothesis is true - values not significantly different from mean/reference
wilcox.test(df.data$SOX2[7:9], mu=mat.mean[4,1])
#null-hypothesis is true - values not significantly different from mean/reference
##NANOG
wilcox.test(df.data$NANOG[1:3], mu=mat.mean[2,2])
#null-hypothesis is true - values not significantly different from mean/reference
wilcox.test(df.data$NANOG[7:9], mu=mat.mean[4,2])
#null-hypothesis is true - values not significantly different from mean/reference
##OCT4
wilcox.test(df.data$OCT4[1:3], mu=mat.mean[2,3])
#null-hypothesis is true - values not significantly different from mean/reference
wilcox.test(df.data$OCT4[7:9], mu=mat.mean[4,3])
#null-hypothesis is true - values not significantly different from mean/reference
