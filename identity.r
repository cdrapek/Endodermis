#####################################################################################
#
# Identity
#
# Calculate cell identity measure (CIM)
#
# Based on "Measuring cell identity in noisy biological systems", Birnbaum and Kussell, 
#
#####################################################################################


RAND_ITER=1000

MIN_USEFUL_SPEC= 0.15

#####################################################################################
#
#  getMarkerList
#
#  returns a list of markers for each tissue
#
#  Arguments
#    ci - spec data structure
#    info - information threshold
#    universe - list of gene names from which to choose markers
#
getMarkerList <- function(ci, info, universe) {
	markers=list()

	ci[[1]][ci[[1]]<0]=0
	for(i in 1:ncol(ci[[1]])) {
		ci_sub = ci[[1]][intersect(rownames(ci[[1]]),universe),i]
		ci_sub = ci_sub[ci_sub> MIN_USEFUL_SPEC]
		ci_max_sub = ci[[2]][names(ci_sub),i]
		cum_ci <- cumsum(ci_sub[order(ci_sub, ci_max_sub, decreasing=TRUE)])
		cum_ci <- cum_ci[1:head(which(cum_ci==max(cum_ci)),n=1)]
		markers[[i]] <- names(which(cum_ci<info))
	}
	
	names(markers)=colnames(ci[[1]])
	markers
}
#####################################################################################
#
#  getRandomBackground
#
#  returns a vector of RAND_ITER length of CIM based on random markers
#   used for significance testing
#
#  Arguments
#    cell - vector of gene expression values
#    ci_for_mark - spec scores
#    universe - list of gene names from which to choose markers
#    marker_num - number of markers to use


getRandomBackground <- function(cell, ci_for_mark, universe, marker_num) {

	rand_id_scores=vector()
	cell = cell[universe]

	for(i in 1:RAND_ITER) {
		# get marker set
		marker_set = sample(1:length(cell), marker_num)
		rand_id_scores[i] = getIdentityScore(cell, ci_for_mark, marker_set)
	}
	rand_id_scores[is.na(rand_id_scores)]=0
	rand_id_scores
}

#####################################################################################
#
#  getIdentityScore
#
#  returns CIMs for a given cell
#
#  Arguments
#    cell - vector of gene expression values
#    ci_for_mark - spec scores for the markers
#    markers - list of markers from getMarkerList

getIdentityScore <- function(cell, ci_for_mark, markers){
	mean(cell[markers]* ci_for_mark) * ((sum(cell[markers]>0))/length(markers))
}
#####################################################################################
#
#  getIdentity
#
#  returns CIMs for a data matrix
#
#  Arguments
#    data - gene expression matrix
#    ci - spec data structure
#    markers - list of markers from getMarkerList
#    returnSig - should significance be calculated
#    universe - what subset of genes should be used for randomizations

getIdentity <- function(data, ci, markers, returnSig=FALSE, universe=c()) {
	
	hs_scoremat <- matrix(nrow=ncol(data), ncol=length(markers))

	colnames(hs_scoremat) <- names(markers)
	rownames(hs_scoremat) <- colnames(data)
	calls <- hs_scoremat
	sig <- calls
	all_markers = unlist(markers)

	for(cell in 1:nrow(hs_scoremat)) {
		markers_cell = data[, cell]
		for(mark in 1:length(markers)) {
			hs_scoremat[cell, mark] = getIdentityScore(data[,cell], ci[[1]][markers[[mark]],mark], markers[[mark]])
			calls[cell, mark] = sum(data[markers[[mark]],cell]>0)
			if(returnSig) {
				sig[cell,mark] <- 1-which(order(c(hs_scoremat[cell,mark], getRandomBackground(markers_cell, ci[[1]][markers[[mark]],mark], universe, length(markers[[mark]]))))==1)/(RAND_ITER+1)
			}
		}
	}

	hs_scoremat_norm <- hs_scoremat
	for(i in 1:nrow(hs_scoremat_norm)) { hs_scoremat_norm[i,] = hs_scoremat_norm[i,]/sum(hs_scoremat_norm[i,]) }
	hs_scoremat_norm[is.nan(hs_scoremat_norm)]=0
	if(returnSig) {
		list(hs_scoremat_norm, hs_scoremat, calls, sig, matrix(nrow=nrow(sig), ncol=ncol(sig), p.adjust(sig, "BH")))
	} else {
		list(hs_scoremat_norm, hs_scoremat)
	}
}

######################################################3

getMaxIDdist <- function(ci, data) {
	ic_mn_list <- list()

	for(i in 1:200) {

		id_tmp = getIdentity(data, ci, getMarkerList(ci, i/2, rownames(data)), returnSig=FALSE)
		ic_mn_list[[i]] = id_tmp[[1]]
		cat(".")
		flush.console()
	}

	maxID = vector()
	for(i in 1:200) {
		maxID[i] = mean(apply(ic_mn_list[[i]],1,max))
	}
	maxID
}

getVariability <- function(ci, data) {
	ic_mn_list <- list()
	for(i in 1:ncol(data)) {
		ic_mn_list[[i]] = matrix(nrow=200, ncol= ncol(ci[[1]]))
		colnames(ic_mn_list[[i]]) = colnames(ci[[1]])
	}

	for(i in 1:200) {
		id_tmp = getIdentity(data, ci, getMarkerList(ci, i/2, rownames(data)), returnSig=FALSE)
		for(x in 1:ncol(data)) {
			ic_mn_list[[x]][i,] = id_tmp[[1]][x,]
		}
		cat(".")
		flush.console()
	}
	for(i in 1:ncol(data)) {
		ic_mn_list[[i]][is.na(ic_mn_list[[i]])]=0
	}

	variability=list()
	for(cell in 1:ncol(data)) {
		variability[[cell]]=vector()
		for(info in 1:198) {
			variability[[cell]][info] = dist(rbind(ic_mn_list[[cell]][info,], ic_mn_list[[cell]][info+2,]))
		}
	}
	variability
}



getIdBasedOnMarkCIM <- function(data, correct, ci,start,end, pval=0.05) {
	correctId <- matrix(nrow = ncol(data), ncol= ncol(ci[[1]]), data=FALSE)
	for(i in 1:length(correct)) {
		correctId[i, correct[i]] = TRUE
	}

	correctRate=vector()
	errorRate=vector()
	markerNum=vector()


	for(i in start:end) {
		m <- getMarkerList(ci, i, rownames(data))
		id <- (0+getIdentity(data, ci, m, TRUE, rownames(data))[[4]]<pval)
		correctRate[(i-start)] = sum(id & correctId)
		errorRate[i-start] = sum ( xor(id, correctId)) 
		markerNum[i-start] = length(unlist(m))/length(m)
		cat(".", correctRate[(i-start)], errorRate[i-start] )
		flush.console()
	}
	list(correctRate, errorRate, markerNum)
}