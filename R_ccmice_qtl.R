
# install DOQTL if not yet installed
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("annotate", "annotationTools", "biomaRt", "Biobase", "corpcor", "GenomicRanges", "hwriter", "MASS", "mclust", "org.Hs.eg.db", "org.Mm.eg.db", "QTLRel", "Rsamtools", "XML"))

# on tower machine, anaconda path disrupted
# export PATH="/usr/local/opt/libxml2/bin:$PATH"
# install.packages("XML")

BiocManager::install("DOQTL", version = "3.8")

dir_ccmice = '~/projects/ccmice';
dir_data = file.path(dir_ccmice, 'data_tower');
# define functions
source(file.path(dir_ccmice, "load_36states.R"))
# call the function
generate_condensed(output.file = file.path(dir_data, "tempCache/founder.probs.B37.Rdata"), input_dir = file.path(dir_data, 'genotype_prob/B37'), temp_dir = file.path(dir_data, 'tempCache/genotypeB37'))
generate_condensed(output.file = file.path(dir_data, "tempCache/founder.probs.B38.Rdata"), input_dir = file.path(dir_data, 'genotype_prob/B38'), temp_dir = file.path(dir_data, 'tempCache/genotypeB38'))


load(file.path(dir_data, "tempCache/founder.probs.B38.Rdata"))


# mm9
# contain all sites of the B38 prob file
load(url('http://csbio.unc.edu/CCstatus/Media/snps.megamuga.Rdata'))
mega_muga = snps
rm(snps)

# mm10
# only contain 63957 sites of the B38 prob
load(url('http://csbio.unc.edu/CCstatus/Media/snps.gigamuga.Rdata'))
giga_mugBa = snps
rm(snps)

temp_marker = read.csv(file.path(dir_data, 'genotype_prob/B37/CC001_Uncb37V01.csv'), header = TRUE)
temp_marker = read.csv(file.path(dir_data, 'genotype_prob/B38/CC001_Uncb38V01.csv'), header = TRUE)
rownames(temp_marker) = temp_marker$marker

dimnames(model.probs)[[3]] = temp_marker$marker
dimnames(model.probs)[[1]] = sapply(strsplit(dimnames(model.probs)[[1]], '_'), '[',  1)


phenotype = read.table(file.path(dir_ccmice, 'data_matlab_tower', 'ccmice_phenotype.txt'), header = TRUE)
phenotype$sex = 'F'
rownames(phenotype) = phenotype$CCStrains
# ccmice_phenotype = as.vector(ccmice_phenotype)


temp_samples = intersect(dimnames(ccmice_phenotype)[[1]], dimnames(model.probs)[[1]])
temp_samples = setdiff(temp_samples, c('CC078', 'CC079', 'CC080', 'CC081', 'CC082', 'CC083'));
# B37 missing CC078-CC083
# B38 missing CC078-CC083 on chrX
temp = apply(model.probs,c(1,3),sum)
temp_sites = apply(temp > 0.99, 2, all)

ccmice_Prob = model.probs[temp_samples,,temp_sites]
ccmice_snps = mega_muga[dimnames(ccmice_Prob)[[3]], c('marker', 'chr', 'pos', 'cM', 'A1', 'A2', 'seq.A', 'seq.B')]
ccmice_snps$chr = temp_marker[dimnames(ccmice_snps)[[1]], 'chromosome']
ccmice_phenotype=phenotype[temp_samples,]
rm(temp)

# export condensed haplotype states
for(i in dimnames(ccmice_Prob)[[1]]){
	tmp = t(ccmice_Prob[i,,])
	tmp = cbind("snp_id"=rownames(tmp), ccmice_snps, tmp)
	write.table(tmp, file = file.path(dir_data, "tempCache/haplotype",  paste(i,"_ccmice_haplotype.tsv",sep="")), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
	}

library('DOQTL')
ccmice_K = kinship.probs(ccmice_Prob)
ccmice_covar = data.frame(sex = as.numeric(ccmice_phenotype$sex == 'M'))
rownames(ccmice_covar) = rownames(ccmice_phenotype)

ccmice_phenotype$EarSwell = scale(ccmice_phenotype$MaximumPCAValue, center = TRUE, scale = TRUE)
ccmice_phenotype$ExpulsionTime = scale(ccmice_phenotype$DateofExpulsion, center = TRUE, scale = TRUE)
ccmice_phenotype$eggcounts_Area= scale(ccmice_phenotype$AUCforeggcounts, center = TRUE, scale = TRUE)
ccmice_phenotype$EarSwell_Area = scale(ccmice_phenotype$AUCforPCA, center = TRUE, scale = TRUE)
qtl = scanone(pheno = ccmice_phenotype, pheno.col = c('EarSwell', 'ExpulsionTime', 'EarSwell_Area', 'eggcounts_Area'), probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps)

# permutation using this one

perm = numeric(0)
# scanone.eqtl uses matrixQTL implementation faster than scanone()
eqtl = scanone.eqtl(ccmice_phenotype[, c('EarSwell', 'EarSwell_Area')], probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps, sex = ccmice_phenotype$sex)
perm = cbind(perm, eqtl)
bit = rownames(ccmice_phenotype)[apply(!is.na(ccmice_phenotype[, c('ExpulsionTime', 'eggcounts_Area')]), 1, all)]
eqtl = scanone.eqtl(ccmice_phenotype[bit, c('ExpulsionTime', 'eggcounts_Area')], probs = ccmice_Prob[bit,,], K = ccmice_K[bit,bit], addcovar = ccmice_covar[bit,, drop = FALSE], snps = ccmice_snps, sex = ccmice_phenotype[bit,, drop = FALSE]$sex)
perm = cbind(perm, eqtl)
nperm = 1000
perm_geno = ccmice_Prob
sample_id = dimnames(ccmice_Prob)[[1]]
for(i in 1:nperm){
	print(paste(i, "of", nperm))
	new.order = sample(1:nrow(ccmice_phenotype))
	dimnames(perm_geno)[[1]] = sample_id[new.order]
	# pheno = (ccmice_phenotype)
	# pheno[,] = (ccmice_phenotype[new.order,,drop = FALSE]) #make sure row names do not permute
	# addcovar = ccmice_covar
	# addcovar[,] = ccmice_covar[new.order,,drop = FALSE]	
	# K = ccmice_K
	# K[,] = ccmice_K[new.order, new.order]
	eqtl = scanone.eqtl(ccmice_phenotype[, c('EarSwell', 'EarSwell_Area')], probs = perm_geno, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps, sex = ccmice_phenotype$sex)
	perm = cbind(perm, eqtl)
	bit = rownames(ccmice_phenotype)[apply(!is.na(ccmice_phenotype[, c('ExpulsionTime', 'eggcounts_Area')]), 1, all)]
	eqtl = scanone.eqtl(ccmice_phenotype[bit, c('ExpulsionTime', 'eggcounts_Area')], probs = perm_geno[bit,,], K = ccmice_K[bit,bit], addcovar = ccmice_covar[bit,, drop = FALSE], snps = ccmice_snps, sex = ccmice_phenotype[bit,, drop = FALSE]$sex)
	perm = cbind(perm, eqtl)
}

perm_max = numeric(0)
for (i in 1:4){
	perm_max = cbind(perm_max, apply(perm[,seq(4+i,ncol(perm),by=4)], 2, max))
}

# qtl_lrs = qtl.LRS(ccmice_phenotype[, 'EarSwell'], probs = ccmice_Prob, addcovar = ccmice_covar, snps = ccmice_snps)
# plot(qtl$EarSwell$lod$A$lod, qtl_lrs$lrs[rownames(qtl$EarSwell$lod$A), 'lod'])

# perms = scanone.perm(pheno = ccmice_phenotype, pheno.col = c('EarSwell', 'ExpulsionTime'), probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps, nperm = 1000)
# thr = quantile(perms, probs = 0.95)

# plot(qtl, sig.thr = thr, main = 'scaled')
# qtl.heatmap(qtl$lod)

# all zero columns, beta is not correct, remove those
qtl_corrected = qtl
for(i in 1:length(qtl)) {
	for(j in 1:length(qtl[[i]]$coef)){
		qtl_corrected[[i]]$coef[[j]][(qtl[[i]]$coef[[j]]) > 5]= NaN
		qtl_corrected[[i]]$coef[[j]][(qtl[[i]]$coef[[j]]) < -5]= NaN
		}
	}

# DOQTL:::plot.doqtl()
source(file.path(dir_ccmice, "html.report_Xin.R"))
html.report_Xin(file.path(dir_ccmice, 'docs', 'QTL'), qtl_corrected[c(1,2)], perms = perm_max[c(1,2),], assoc = FALSE)
html.report_Xin(file.path(dir_ccmice, 'docs', 'QTL'), qtl_corrected[c(3,4)], perms = perm_max[c(3,4),], assoc = FALSE)

save.image(file=file.path(dir_data, "tempCache", "ccmice_10202019.RData"))
savehistory("~/mac_hdd/ccmice/tempCache/ccmice_apr16.Rhistory")







