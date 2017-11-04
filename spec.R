#####################################################################################
#
# Spec
#
# Calculate information scores
#
# Based on "Measuring cell identity in noisy biological systems", Birnbaum and Kussell, 
#
# Spec score is calculate as followed:
#
# I(level) = 1 + sum(P(type|level)logn(P(type|level)
# Spec(type) = sum(I(level)P(level|type))
#

#####################################################################################
# specy
#
# Calculate spec score for a given gene
# 
# Arguments
#    gene - vector of gene expression vector for samples
#    header - vector of sample codes
#    binsize - size of the bin
# Returns
#    list of Spec scores gefor the gene, and the mean expression/binsize for each score

specy <- function(gene, header, binsize) {
    ntypes <- length(unique(as.vector(header)))

    bylist=list()
    bylist[[1]]=header
    meanValues_t = aggregate(gene, bylist,median)
    meanValues= meanValues_t[,2]
    names(meanValues)=meanValues_t[,1]

    bins <- cut(gene, breaks=c(min(gene)-1, binsize), labels=F)
    bins[is.na(bins)] <- 2 # For values above hmax

    tab <- as.matrix(table(as.data.frame(cbind(header, bins))))

    # If not all bins are represented, fill a 0 matrix instead
    if (ncol(tab) < 2) {
        binned <- matrix(0,ntypes,2)
        for (i in colnames(tab)) {
            binned[,as.numeric(i)] <- tab[,i]
        }
    } else {
        binned <- tab
    }

    pleveltype <- binned/rowSums(binned)
    pleveltype[is.na(pleveltype)] <- 0
    ptypelevel <- t(pleveltype)/colSums(pleveltype)
    ptypelevel[is.na(ptypelevel)] <- 0
    logNptypelevel <- log(ptypelevel, ntypes)
    logNptypelevel[is.infinite(logNptypelevel)] = 0 # For log0 values
    #logNptypelevel[is.na(logNptypelevel)] = 0 # For log0 values
    ilevel <- 1 + rowSums(ptypelevel*logNptypelevel)
    ilevel[is.na(ilevel)] <- 0
    
    specy <- pleveltype%*%ilevel
    # Information is in the negative if the info is in the bottom bins
 #   infosign = (as.integer(aggregate (bins, by=list(header), function(x) { mean(x)})[,2] < (1+2)/2) -1)* (- 2) -2+1
    infosign = as.integer(aggregate (gene, by=list(header), function(x) { mean(x)})[,2] > mean(gene)) -1 
    # we do not use "absent" information for the time being
    specy[infosign<0]=0
    #
    return (list(specy, meanValues/binsize))

}


#######################################################################3
#
#  optimizebinsize
#
#  Identify the optimal binsize for a given gene
#
#  Arguments
#    gene - vector of gene expression vector for samples
#    header - vector of sample codes
#    cuts - number of bins (0 to not use binning) (l)
#    distshape - max bin to be an acceptable background bin (u)
#  Returns:
#	optimal binsize

optimizebinsize <- function(gene, header, cuts=0, distshape=0) {

    background_cutoff=0

    if(cuts) {
	    # gene should conform to a power-law dist
	    histo <- (table(cut(gene, breaks=seq(min(gene),max(gene),(max(gene)-min(gene))/cuts)))) >= ceiling(length(header)/cuts)
	    if(length(which(histo))>0) {
		    background_cutoff = max(which(histo))
	    } else {
		background_cutoff=length(histo)
	    }
    }

    # If gene does not conform, abort
    if(background_cutoff <= distshape ) {
	ntypes <- length(unique(as.vector(header)))
	minbin = min(gene) + (max(gene)-min(gene))*0.1
	maxbin = max(gene) - (max(gene)-min(gene))*0.1
	if(cuts) {# if using distribution shape filtering, just go for the cuts
		minbin = min(gene)+(max(gene)-min(gene))/cuts*(background_cutoff-1.5)
		maxbin = min(gene)+(max(gene)-min(gene))/cuts*(background_cutoff+0.5)
	}

	binsizes <- seq(minbin, maxbin, (maxbin-minbin)/50)
	infosum = matrix(nrow=length(binsizes), ncol=ntypes, data=0)
	rownames(infosum) = binsizes
	for(i in 1:length(binsizes)) {
		infosum[i,]= specy(gene, header, binsizes[i])[[1]][,1]
	}
	# return maximizie information for the highest sample
	binsizes[ match(max(infosum),apply(infosum,1,max))]
    } else {
	NA
    }
}

#######################################################################
#
# getAllSpec
#
# accepts a data table with the first row (Group) marks the group codes
#
# Arguments
#   data - expression matrix
#   header - vector of sample codes
#   medianfilter - filter gene whose median is above this
#    cuts - number of bins (0 to not use binning) (l)
#    distshape - max bin to be an acceptable background bin (u)
#
# Returns
#    The Spec data structure
#

getAllSpec <- function(data, header, medianfilter=0, cuts= FALSE, distshape=0) {

	# Settings
	float_prec <- 5
	if(distshape==0) {
		distshape = ceiling(cuts*0.33)
	}
	if(medianfilter==0) {
		medianfilter=max(data)
	}

	# Get Spec values row by row (gene by gene)
	specs <- matrix(nrow=nrow(data), ncol= length(unique(as.vector(header))), data=0)
	specMeanValues <- matrix(nrow=nrow(data), ncol= length(unique(as.vector(header))), data=0)

	for (gene in 1:nrow(data)) {
	    binsize <- optimizebinsize(data[gene,],header, cuts, distshape)
	    if(!is.na(binsize) && median(data[gene,])<medianfilter) {
		specg <- specy(data[gene,], header, binsize)
		specs[gene,] = specg[[1]]
		specMeanValues[gene,] = specg[[2]]
	    } else {
		specs[gene,] = rep(0, ncol(specs))
		specMeanValues[gene,] = rep(0, ncol(specs))
	    }

	    if( (gene %% 10)==0) {
		    cat(".")
		    flush.console()
	    }
	    if( (gene %% 1000)==0) {
		cat(gene/nrow(data)*100, "%  \n")
	    }
	}

	rownames(specs) <- rownames(data)
	colnames(specs) <- unique(as.vector(header))
	rownames(specMeanValues) <- rownames(data)
	colnames(specMeanValues) <- unique(as.vector(header))

	list(round(specs, float_prec), specMeanValues)
}
