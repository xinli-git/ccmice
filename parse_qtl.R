

parse_qtl = function( qtl, ccmice_hap ){
qtl_table = NULL
for( i in names(qtl) ) {
	EarSwell_qtl = rbind(qtl[[i]]$lod$A, qtl[[i]]$lod$X);
	pvalue = EarSwell_qtl[,5:9]
	dimnames(pvalue)[[2]] = paste(dimnames(pvalue)[[2]], i, sep="_")
	if( is.null(qtl_table) ) {
		qtl_table = EarSwell_qtl[,1:4]
	}
	qtl_table = cbind(qtl_table, pvalue[dimnames(qtl_table)[[1]],])
}
qtl_table = cbind(qtl_table, ccmice_hap[dimnames(qtl_table)[[1]],])
return(qtl_table)
}


