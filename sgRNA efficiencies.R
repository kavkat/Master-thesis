rm(list=ls())
library(extrafont)
library(showtext)
font_add_google("Roboto Slab", "roboto_slab")
showtext_auto()

#AxD1 pX459 

#Indels
par(mfrow=c(2,1), family="roboto_slab")
indels1.59 <- read.csv('AxD1.59_sgRNA.csv')
indels1.59$sgRNA <- factor(indels1.59$sgRNA)
boxplot(indels1.59$Indels~indels1.59$sgRNA, xlab='sgRNA', ylab=expression(paste(Delta)), col='wheat', main='\nComparision of Insertion/Deletion lengths')

#Efficiency
sgRNA1.59 <- c(14, 13, 12, 11)
eff1.59 <- c(83.33, 80, 50, 60)
barplot(eff1.59~sgRNA1.59, col='lightblue', ylab= 'sgRNA', xlab='Efficiency(%)', xlim=c(0,100), horiz=TRUE, main='DSB induction efficiency')

mtext("sgRNA efficiency analysis for CRISPR-Cas9 gene editing on LU AxD1", side=3, line=-1.5, font=2, cex=1.5, outer=TRUE)

#AxD1 pX462 

#Indels
par(mfrow=c(2,1), family="roboto_slab")
indels1.62 <- read.csv('AxD1.62_sgRNA.csv')
indels1.62$sgRNA <- factor(indels1.62$sgRNA)
boxplot(indels1.62$Indels~indels1.62$sgRNA, xlab='sgRNA combinations', ylab=expression(paste(Delta)), col='tan', main='\nComparision of Insertion/Deletion lengths')

#Efficiency
sgRNA1.62 <- c('13:14', '13:12', '11:14', '11:12')
eff1.62 <- c(75, 100, 88.89, 87.5)
barplot(eff1.62~sgRNA1.62, col='seagreen3', ylab= 'sgRNA combinations', xlab='Efficiency(%)', xlim=c(0,100), horiz=TRUE, main='DSB induction efficiency')

mtext("sgRNA efficiency analysis for CRISPR-nCas9 gene editing on LU AxD1", side=3, line=-1.5, font=2, cex=1.5, outer=TRUE)

#AxD9 pX459
par(mfrow=c(2,1), family="roboto_slab")
indels9.59 <- read.csv('AxD9.59_sgRNA.csv')
indels9.59$sgRNA <- factor(indels9.59$sgRNA)
boxplot(indels9.59$Indels~indels9.59$sgRNA, xlab='sgRNA', ylab=expression(paste(Delta)), col='peru', main='\nComparision of Insertion/Deletion lengths')

#Efficiency
sgRNA9.59 <- c(94, 93, 92, 91)
eff9.59 <- c(88.89, 55.56, 0, 100)
barplot(eff9.59~sgRNA9.59, col='steelblue', ylab= 'sgRNA', xlab='Efficiency(%)', xlim=c(0,100), horiz=TRUE, main='DSB induction efficiency')

mtext("sgRNA efficiency analysis for CRISPR-Cas9 gene editing on LU AxD9", side=3, line=-1.5, font=2, cex=1.5, outer=TRUE)

#AxD9 pX462 

#Indels
par(mfrow=c(2,1), family="roboto_slab")
indels9.62 <- read.csv('AxD9.62_sgRNA.csv')
indels9.62$sgRNA <- factor(indels9.62$sgRNA)
boxplot(indels9.62$Indels~indels9.62$sgRNA, xlab='sgRNA combinations', ylab=expression(paste(Delta)), col='sandybrown', main='\nComparision of Insertion/Deletion lengths')

#Efficiency
sgRNA9.62 <- c('93:94', '93:92', '91:94', '91:92')
eff9.62 <- c(62.5, 20, 100, 14.29)
barplot(eff9.62~sgRNA9.62, col='paleturquoise3', ylab= 'sgRNA combinations', xlab='Efficiency(%)', xlim=c(0,100), horiz=TRUE, main='DSB induction efficiency')

mtext("sgRNA efficiency analysis for CRISPR-nCas9 gene editing on LU AxD9", side=3, line=-1.5, font=2, cex=1.5, outer=TRUE)
