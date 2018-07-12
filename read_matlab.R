
# source("/Users/xinli/Dropbox/R/ccmice/R_ccmice.R")

ccmice_marker = read.table('/Volumes/Mac HDD/matlab/history/ccmice/ccmice_marker.txt', header = TRUE)
rownames(ccmice_marker) = ccmice_marker$marker 

ccmice_allele = read.table('/Volumes/Mac HDD/matlab/history/ccmice/ccmice_allele.txt', header = TRUE)

ccmice_sample = read.table('/Volumes/Mac HDD/matlab/history/ccmice/ccmice_sample.txt', header = TRUE)

ccmice_phenotype = read.table('/Volumes/Mac HDD/matlab/history/ccmice/ccmice_phenotype.txt', header = TRUE)
ccmice_phenotype$sex = 'F'
rownames(ccmice_phenotype) = ccmice_phenotype$CCStrains
# ccmice_phenotype = as.vector(ccmice_phenotype)

load(url('http://csbio.unc.edu/CCstatus/Media/snps.megamuga.Rdata'))
mega_muga = snps
rm(snps)

sum(ccmice_marker$marker %in% mega_muga$marker)

library('R.matlab')
ccmice_haplotype = readMat('/Volumes/Mac HDD/matlab/history/ccmice/ccmice_matrix_04132018.mat', fixNames = FALSE, verbose = TRUE)
ccmice_haplotype = ccmice_haplotype$ccmice_matrix[[4]]
dimnames(ccmice_haplotype)[[1]] = ccmice_marker$marker
dimnames(ccmice_haplotype)[[2]] = ccmice_allele$alleles
dimnames(ccmice_haplotype)[[3]] = ccmice_sample$sample_id
ccmice_haplotype = aperm(ccmice_haplotype,perm = c(3,2,1))

# K = kinship.probs(ccmice_haplotype[,,seq(1,dim(ccmice_haplotype)[3],10)])
temp = apply(ccmice_haplotype,c(1,3),sum)
temp = apply(temp > 0.99, 2, all)
ccmice_Prob = ccmice_haplotype[dimnames(ccmice_phenotype)[[1]],,temp]
ccmice_snps = mega_muga[dimnames(ccmice_Prob)[[3]], c('marker', 'chr', 'pos', 'cM', 'A1', 'A2', 'seq.A', 'seq.B')]
ccmice_snps$chr = ccmice_marker[dimnames(ccmice_snps)[[1]], 'chromosome']
rm(temp)


library('DOQTL')
ccmice_K = kinship.probs(ccmice_Prob)
ccmice_covar = data.frame(sex = as.numeric(ccmice_phenotype$sex == 'M'))
rownames(ccmice_covar) = rownames(ccmice_phenotype)

ccmice_phenotype$EarSwell = scale(ccmice_phenotype$MaximumValue, center = TRUE, scale = TRUE)
ccmice_phenotype$ExpulsionTime = scale(ccmice_phenotype$DateofExpulsion, center = TRUE, scale = TRUE)
qtl = scanone(pheno = ccmice_phenotype, pheno.col = c('EarSwell', 'ExpulsionTime'), probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps)

# permutation using this one
scanone.eqtl(ccmice_phenotype[, c('EarSwell', 'ExpulsionTime')], probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps, sex = ccmice_phenotype$sex)


perms = scanone.perm(pheno = ccmice_phenotype, pheno.col = c('EarSwell', 'ExpulsionTime'), probs = ccmice_Prob, K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps, nperm = 1000)
thr = quantile(perms, probs = 0.95)
plot(qtl, sig.thr = thr, main = 'scaled')
qtl.heatmap(qtl$lod)


# remove rank -deficient SNPs
bit1 = (abs(qtl$coef$A[,3:9]) < 5) & !is.na(qtl$coef$A[,3:9])  
bit2 = apply(bit1, 1, all)
bitA = rownames(qtl$coef$A[bit2,])


bit1 = (abs(qtl$coef$X[,3:9]) < 5) & !is.na(qtl$coef$X[,3:9])   
bit2 = apply(bit1, 1, all)
bitX = rownames(qtl$coef$X[bit2,])
bit = c(bitA,bitX)
qtl2 = scanone(pheno = ccmice_phenotype, pheno.col = 'scaled', probs = ccmice_Prob[,,bit], K = ccmice_K, addcovar = ccmice_covar, snps = ccmice_snps[bit,])

% create.html.page('/Volumes/Mac HDD/ccmice/QTL/', qtl, 'scaled', perms, assoc = FALSE)
% html.report('/Volumes/Mac HDD/ccmice/QTL/', qtl, perms, assoc = FALSE)


pdf('/Volumes/Mac HDD/ccmice/QTL/ccmice_allChr.pdf')
plot(qtl, sig.thr = thr, main = 'scaled')
for(i in c(1:19, 'X'))
{
	coefplot(qtl2, chr = i, sex = 'F')
}
dev.off()






