# **************************************************************** #
# snp annotation for the suggestive SNPs
# compare with previous studies
#---------------------------------------------------------------- #

setwd('hgb')

counter = 1
for (chr in 1:22) {
	info <- read.csv(paste('imputed_info/info_files_filtered/hybrid_info', chr, '_filtered.csv', sep = ''), as.is = T)
	results <- read.csv(paste('gwa_results_pval_cutoff/chr', chr, '_gwa_out.csv', sep = ''), as.is = T)
	
	results_info <- merge(results, info, by.x = 'snp', by.y = 'markname')

	# cleaning
	results_info$MAF <- as.numeric(results_info$MAF)
	
	if (counter == 1) {
		all_info <- results_info
	} else {
		all_info <- rbind(all_info, results_info)		
	}
	counter = counter + 1
}

dim(all_info)[1] 
head(all_info[order(all_info$Pval), ])

# plot the results
plot <- all_info[,c('snp', 'Pval', 'Beta', 'SE', 'chrom', 'pos', 'MAF')]
colnames(plot) <- c('SNP', 'P', 'Beta', 'SE', 'CHR', 'BP', 'MAF')

source('common_files/manhattan_GGD.R')
manhattan(plot, main = 'HGB', title = 'HGB', colors = c("black","#666666","#CC6600"))
manhattan(plot[plot$MAF >= 0.01, ], main = 'HGB_maf1', title = 'HGB-0.01', colors = c("black","#666666","#CC6600"))
manhattan(plot[plot$MAF >= 0.05, ], main = 'HGB_maf5', title = 'HGB-0.05',colors = c("black","#666666","#CC6600"))

# Suggestive SNPs; Pvalue < 5e-6
sug <- all_info[all_info$Pval < 5e-6, ]
dim(sug)[1] 
rm(list = c('qq', 'plot', 'all_info', 'results', 'results_info', 'chr', 'counter', 'info', 'manhattan'))

# change beta sign if required
table(sug$BetaSwitchSign)
sug$Beta <- ifelse((sug$BetaSwitchSign == 'No'), sug$Beta, -(sug$Beta)) # if 'No' then dont' change beta sign; for genotyped SNPs it is already fine
sug$logp <- -log10(sug$Pval) # calculate logp value

# ----------------------------- Get the cyto band --------------------------------------------------------- #
cyto <- read.table('annovar/humandb/hg19_cytoBand.txt', sep = '\t', as.is = T)
for (i in 1:dim(sug)[1]) sug$Region[i] <- cyto$V4[cyto$V1 == paste('chr', sug$chrom[i], sep = '') & sug$pos[i] > cyto$V2 & sug$pos[i] < cyto$V3]

sug$Region <- paste(sug$chrom, sug$Region, sep = '')

# +++++++++++++++++++ Add annotation for SNPs +++++++++++++++++++ #
#	within gene; 60kb +/- gene
# 	gene information: HUGO identifiers (17,000 genes) hg19 positions
# ---------------- hugo gene annotations ---------------------- #
hugo <- read.table('common_files/hugo_genes', as.is = T)
colnames(hugo) <- c('chr', 'start', 'end', 'gene')
head(hugo)
# ---------------------- gene annotation function ------------- #
# window: for within gene have window = 0 or can modify it (60kb either sides)
geneAnnotation <- function(df, window) {
	df$RefGene <- NA
	for (row in 1:dim(df)[1]) {
		tmp <- hugo[hugo$chr == df$chrom[row],]
		pos <- df$pos[row]
		tmp$start <- tmp$start - window
		tmp$end <- tmp$end + window
		gene <- c()
		for (i in 1:dim(tmp)[1]) {
			if (tmp$start[i] < pos) {
				if(tmp$end[i] > pos) {
					gene <- c(gene, tmp$gene[i])
				}
			}
		}
		gene <- paste(gene, collapse = ',')
		df$RefGene[row] <- gene
	}
	return(df)
}

# ------------------------- Adding Gene Annotations------------------------------------ #
sug <- geneAnnotation(sug, 60000)
colnames(sug)[dim(sug)[2]] <- 'RefGenes60'
sug <- geneAnnotation(sug, 0)
for (i in 1:dim(sug)[1]) sug$RefGenes60[i] <- paste(unique(unlist(strsplit(c(paste(sug$RefGenes60[i], collapse = ',')), split = ','))), collapse = '; ')
for (i in 1:dim(sug)[1]) sug$RefGene[i] <- paste(unique(unlist(strsplit(c(paste(sug$RefGene[i], collapse = ',')), split = ','))), collapse = '; ')
# -------------------------------------------------------------------------------------- #

sug <- sug[,c('snp', 'Region', 'pos', 'minor_major_allele', 'MAF', 'Beta', 'SE', 'Pval', 'imputed', 'logp', 'chrom', 'RefGenes60', 'RefGene')]

sug <- sug[order(sug$chrom, sug$pos), ]
sig <- sug[sug$Pval < 5e-8,]

write.csv(sig, file = 'sig_info.csv', row.names = F, quote = F)
write.csv(sug, file = 'sug_info.csv', row.names = F, quote = F)

# ------------------------- Previous Studies ------------------------------------------ #
sig <- read.csv('sig_info.csv', as.is = T)

pre <- read.delim('gwa_catalog/nhgri_ebi/previous_hits.txt', as.is = T, header = T)

PreviousStudies <- function(df, window) {
	pre_df <- pre[pre$Pre_CHR == 25, ] # make empty data set
	pre_df$sugSNP <- as.character(); pre_df$dist <- as.numeric(); pre_df$sugSNP_Pval <- as.numeric(); pre_df$sugSNP_MAF <- as.numeric()
	for (row in 1:dim(df)[1]) {
		pos <- df$pos[row]
		tmp <- pre[pre$Pre_CHR == df$chrom[row] & pre$grch37_pos > (pos - window) & pre$grch37_pos < (pos + window),]
		if (nrow(tmp) != 0) {
			tmp$sugSNP <- gsub('X', '', df$snp[row])
			tmp$dist <- round((pos - tmp$grch37_pos)/1000, digits = 2)
			tmp$sugSNP_Pval <- df$Pval[row]		
			tmp$sugSNP_MAF <- df$MAF[row]
			pre_df <- rbind(pre_df, tmp)
		}
	}
	return(pre_df)
}
pre_hits <- PreviousStudies(sig, window = 10000)
dim(pre_hits)

write.csv(pre_hits, file = 'pre_hits_sug_0.01.csv', quote = F, row.names = F, na = '')


##########################################################################################################
