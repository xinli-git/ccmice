
# test other phenotypes

library(stringr)

ccmice_meta = read.table('~/mac_hdd/ccmice/metaData/allCClines.txt', header = TRUE, sep = "\t")
ccmice_meta$sample_id = str_match(ccmice_meta$name, "(CC\\d+)\\/")[,2]
rownames(ccmice_meta) = ccmice_meta$sample_id


samples = intersect(dimnames(model.probs)[[1]], dimnames(ccmice_meta)[[1]])


ccmice_phenotype = ccmice_meta[samples, c('headSpot', 'avg_LitterSize')]
ccmice_phenotype$white_spot = ccmice_phenotype$headSpot == 'yes'


ccmice_phenotype$sex = 'F'


ccmice_phenotype$white_spot = scale(ccmice_phenotype$white_spot, center = TRUE, scale = TRUE)
ccmice_phenotype$avg_LitterSize = scale(ccmice_phenotype$avg_LitterSize, center = TRUE, scale = TRUE)
qtl = scanone(pheno = ccmice_phenotype, pheno.col = c('white_spot', 'avg_LitterSize'), probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps)

perm = numeric(0)
bit = rownames(ccmice_phenotype)[apply(!is.na(ccmice_phenotype[, c('white_spot', 'avg_LitterSize')]), 1, all)]
eqtl = scanone.eqtl(ccmice_phenotype[bit, c('white_spot', 'avg_LitterSize')], probs = ccmice_Prob[bit,,], K = ccmice_K[bit,bit], addcovar = ccmice_covar[bit,, drop = FALSE], snps = ccmice_snps, sex = ccmice_phenotype[bit,, drop = FALSE]$sex)
perm = cbind(perm, eqtl)
nperm = 100
perm_geno = ccmice_Prob
sample_id = dimnames(ccmice_Prob)[[1]]
for(i in 1:nperm){
	print(paste(i, "of", nperm))
	new.order = sample(1:nrow(ccmice_phenotype))
	dimnames(perm_geno)[[1]] = sample_id[new.order]
	bit = rownames(ccmice_phenotype)[apply(!is.na(ccmice_phenotype[, c('white_spot', 'avg_LitterSize')]), 1, all)]
	eqtl = scanone.eqtl(ccmice_phenotype[bit, c('white_spot', 'avg_LitterSize')], probs = perm_geno[bit,,], K = ccmice_K[bit,bit], addcovar = ccmice_covar[bit,, drop = FALSE], snps = ccmice_snps, sex = ccmice_phenotype[bit,, drop = FALSE]$sex)
	perm = cbind(perm, eqtl)
}

perm_max = numeric(0)
for (i in 1:2){
	perm_max = cbind(perm_max, apply(perm[,seq(2+i,ncol(perm),by=2)], 2, max))
}

source("~/Dropbox/R/ccmice/html.report_Xin.R")
html.report_Xin('~/mac_hdd/ccmice/QTL/', qtl[c(1,2)], perms = perm_max[c(1,2),], assoc = FALSE)


