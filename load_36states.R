


# load data from 36 state genotype probability
# condense it to additive 8 states

# 36 states probability downloaded from 
# http://csbio.unc.edu/CCstatus/index.py?run=FounderProbs
# resolved genotypes at http://csbio.unc.edu/CCstatus/CCGenomes/
# http://csbio.unc.edu/CCstatus/CCGenomes/#genotypes

library('DOQTL')

create.Rdata.files_36states = function(prob.files, cross = "DO") {

	samples = basename(prob.files)
	samples = sub("\\.csv$", "", samples)
	print(samples)

  for(i in 1:length(prob.files)) {

    print(prob.files[i])
    prsmth = read.csv(prob.files[i])
   	print(prsmth[1:10,])
    prsmth = prsmth[,4:ncol(prsmth)]
    # prsmth = exp(as.matrix(prsmth))
    prsmth = (as.matrix(prsmth))
    class(prsmth) = c("genoprobs", class(prsmth))
    attr(prsmth, "cross") = cross
    save(prsmth, file = paste('~/mac_hdd/ccmice/tempCache/genotype', samples[i], '.genotype.probs.Rdata', sep = "") )
    print(file)

  } # for(i)

}

generate_condensed = function(output.file = "~/mac_hdd/ccmice/tempCache/founder.probs.Rdata"){

	files = dir('~/mac_hdd/ccmice/genotype_prob/B37', pattern = ".csv$", full.names = TRUE)
	create.Rdata.files_36states(files, cross = "CC")
	condense.model.probs(path = "~/mac_hdd/ccmice/tempCache/genotype", write = output.file, model = c("additive"), cross = "CC")

}


# R format, functions defined before calling
# source("~/Dropbox/projects/ccmice/load_36states.R")

