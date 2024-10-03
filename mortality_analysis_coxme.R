# ------------------------------------------------------------------ #
# Mortality analysis using coxme
# Input files: 
# lng: Endophenotype file: all_variables_with_PCs.csv
#		Follow-up: followup_death.csv
#		Generation Info: triplet.csv
# ------------------------------------------------------------------- #
setwd('~/Desktop/FILE')
library(coxme) # load library
library(survival)
library(kinship2)
library(AUC)
# ----------------------------------------------------------------------------------------------------- #
# Load the endophenotype data and sub-setting it
data <- read.csv('', as.is = T)
pcs <- subset(data, select = c(subject, X_AGE, sex, PC1, PC2, PC3, PC4, PC5))

# mortality data
followup <- read.csv('followup_death_cleaned.csv', as.is = T) 
tmp <- merge(pcs, followup[ ,1:7], by = 'subject')

# get the generation info
triplet <- read.csv('triplet.csv', as.is = T)
tmp <- merge(tmp, triplet[,c('subject', 'pedid', 'dadsubj', 'momsubj', 'control', 'gen', 'relative')], by = 'subject')

# ------------------ get number of deaths and follow-up time ---------------- #
dim(tmp)
ndeath <- merge(pcs, triplet[,c('subject', 'pedid', 'dadsubj', 'momsubj', 'control', 'gen', 'relative')], by = 'subject')
ndeath <- merge(ndeath, followup, by = 'subject')
all <- ndeath 
proband <- ndeath[ndeath$gen == 2 & ndeath$relative == 1, ] 
offspring <- ndeath[ndeath$gen == 3 & ndeath$control == 0, ] 

# ---------------------------------------------------------------------------- #

# calculate kinship matrix
kmat <- makekinship(tmp$pedid, tmp$subject, tmp$dadsubj, tmp$momsubj)

# run survival analysis for probands using coxme
lng_fit0 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit1 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC1 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit2 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC2 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit3 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC3 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit4 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC4 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit5 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC5 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))
lng_fit6 <- coxme(Surv(days_since_enrollment, death) ~ X_AGE + sex + PC1 + PC2 + PC3 + PC4 + PC5 + (1 | subject), data=tmp, varlist = coxmeMlist(2*kmat, rescale = F), subset = (gen == 2))

lng_results <- data.frame(cof = as.numeric(), se = as.numeric(), cofExp = as.numeric(), lower = as.numeric(), upper = as.numeric(), pvalue = as.numeric(), rowNames = as.character())
nullrow <-  data.frame(cof = NA, se = NA, cofExp = NA, lower = NA, upper = NA, pvalue = NA, rowNames = NA)

dataList <- list(lng_fit0 = lng_fit0, lng_fit1 = lng_fit1, lng_fit2 = lng_fit2, lng_fit3 = lng_fit3, lng_fit4 = lng_fit4, lng_fit5 = lng_fit5, lng_fit6 = lng_fit6)

for (i in 1:length(dataList)) {
	df <- data.frame(cof = dataList[[i]]$coefficients, se = sqrt(diag(vcov(dataList[[i]]))))
	df$cofExp <- round(exp(df$cof), digits = 4)
	df$lower <- round(exp(df$cof - (1.96*df$se)), digits = 4)
	df$upper <- round(exp(df$cof + (1.96*df$se)), digits = 4)
	df$pvalue <- 2*pnorm(-abs(df$cof/df$se))
	df$rowNames <- rownames(df)
	lng_results <- rbind(lng_results, df)
	lng_results <- rbind(lng_results, nullrow)
}

lng_results$cof <- round(lng_results$cof, digits = 4)
lng_results$se <- round(lng_results$se, digits = 4)
head(lng_results)
write.csv(lng_results[,c(7, 1:6)], file = 'lng_coxme_results.csv', quote = F, row.names = F, na = '')

# Calculating attenuation of age HR
age_attenu <- data.frame(fit = names(dataList), Att = NA)
age_attenu$Att[1] <- 'Ref'
for (i in 2:length(dataList)) age_attenu$Att[i] <- (round((((1-exp(lng_fit0$coefficient[1])) - (1-exp(dataList[[i]]$coefficient[1])))/(1-exp(lng_fit0$coefficient[1]))*100), digit = 1))
age_attenu



