

parse_qtl = function( qtl, ccmice_hap ){
qtl_table = NULL
for( i in names(qtl) ) {
	EarSwell_lod = rbind(qtl[[i]]$lod$A, qtl[[i]]$lod$X);
	EarSwell_coef = rbind(qtl[[i]]$coef$A, qtl[[i]]$coef$X);
	pvalue = EarSwell_lod[,5:9]
	coef = EarSwell_coef
	dimnames(pvalue)[[2]] = paste(dimnames(pvalue)[[2]], i, sep="_")
	dimnames(coef)[[2]] = paste(dimnames(coef)[[2]], i, sep="_")
	if( is.null(qtl_table) ) {
		qtl_table = EarSwell_lod[,1:4]
	}
	qtl_table = cbind(qtl_table, pvalue[dimnames(qtl_table)[[1]],], coef[dimnames(coef)[[1]],])
}
qtl_table = cbind(qtl_table, ccmice_hap[dimnames(qtl_table)[[1]],])
return(qtl_table)
}


