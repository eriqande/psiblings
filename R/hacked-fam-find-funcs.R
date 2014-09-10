
#' @include "pseudo.R"
NULL


# this just returns a matrix of pairs that have edges
GetAdjacencies <- function(fit, c.tuning=1) {
	n = fit$n					# Sample size
  N = n*(n-1)/2					# Number of pairwise comparisons
	Mpost = fit$posterior				# Matrix of posterior probabilities
	# Selects the top p_fs-fraction of list sorted according to P_posterior(full sib) 
  #kk=Mpost[Mpost[,4] > quantile(Mpost[,4],1.0-c.tuning[1]*p_fs),]
  #	n1 = nrow(kk)
  n1=sum(Mpost[,4])
  kk=Mpost[Mpost[,4] > quantile(Mpost[,4],1.0-c.tuning[1]*n1/N),]
  kk[,c(1,2)]
}


# give it a file name and it spits back a list of inferred full sibling groups
GetFullSibs <- function(File) {
	tmp = as.matrix(read.table(File,row=1))
	#tmp = as.matrix(read.table("simulations_eric/allhalf_r001_d02m01.txt",row=1))
	S = ncol(tmp)/2

	A_eric= tmp[,2*(1:S)-1]
	B_eric= tmp[,2*(1:S)]

	fit = pseudo(A=A_eric,B=B_eric,relations=c("full-sibs","half-sibs"),eps_error=0.0)

	adj = GetAdjacencies(fit)

	# get the full sibling output
	fs = full_siblings(fit)

	list(fs=fs, ids=row.names(tmp), n=nrow(tmp), edges=adj, fit=fit)
}


# take the output (fs) and input (tmp) to pseudo2 and make a data frame of the
# results.  Also make an "relationship matrix" that has a 1 for every pair that
# is deemed unrelated and a 12 for everyone that is deemed a full sibling.  
# this can get passed to family finder's GH cut tree algorithm.
fs2df <- function(fs, ids, n, numloc=NA, code=NA, opt.str=NA) {
	if(is.null(fs)) {
		all.idx <- 1:n
	} else {
		fs <- fs[names(rev(sort(sapply(fs,length))))] # sort the sibling groups on the basis of size
		singletons<-setdiff(ids, unlist(fs))
		idxs <- 1:n
		names(idxs) <- ids
		fs.idx <- lapply(fs,function(x) idxs[x])
		
		# before putting the singletons on there, we subscript a matrix with the pairs in inferred sibships of size
		# greater than 1
		RM <- matrix(1, nrow=n, ncol=n) # this initializes it to 1's
		for(x in fs.idx) {
			tt <- expand.grid(x,x)  # all pairs make by the indices of those in the sibgroup
  		tt <- as.matrix(tt[tt[,1] != tt[,2], ])   # drop the pairs that are (1,1) or (2,2) or (m,m) for any m, and make it a matrix
  		RM[tt] <- 12  # change those elements of RM corresponding to those pairs
		}
		
		
		if(length(singletons)) {
			sing.idx <- as.list(sort(idxs[singletons]))
			all.idx <- c(fs.idx,sing.idx)
		} else {
			all.idx <- fs.idx
		}
		
  }
  
  # down here we make this relationship matrix, which is going to be passed to
  # family finder as "results"
#  RM <- matrix(1, nrow=n, ncol=n) # this initializes it to 1's
	# this puts the 12's in where they should be based on all the sibships
 # junk <- lapply(fs.idx, function(x) {
 # 	if(length(x)>1) {
 # 			tt <- expand.grid(x,x)  # all pairs make by the indices of those in the sibgroup
 # 			tt <- as.matrix(tt[tt[,1] != tt[,2], ])   # drop the pairs that are (1,1) or (2,2) or (m,m) for any m, and make it a matrix
 # 			RM[tt] <- 12  # change those elements of RM corresponding to those pairs
 # 		}
 # 	}
 # )
  
	SibSize <- sapply(all.idx, length); names(SibSize)<-NULL
	SibSet <- sapply(all.idx, function(x) paste(x, collapse=",")); names(SibSet)<-NULL
	
	data.frame(code, numloc, opt.str, SibSize, SibSet)
}



## here is a function that takes the list that comes out of 
## the full_siblings function and sends it to my hacked version
## of family finder and sends back data like fs2df
fs2df.ff <- function(y, numloc=NA, code=NA, opt.str=NA)  {
	# the name of the function in the load table is horribly mangled!!
	FF <-  .C("famfind", 
									as.integer(y$n), 
									as.integer(y$edges[,1]), 
									as.integer(y$edges[,2]), 
									as.integer(nrow(y$edges)),
									fam.vec=integer(y$n),
									fam.start=integer(y$n),
									fam.end=integer(y$n),
									num.fam=integer(1)
									)
	# and now make that into a list, and turn it into strings
	ff.list <- lapply(1:FF$num.fam, function(x) FF$fam.vec[FF$fam.start[x]:FF$fam.end[x]])
	ff.list <- ff.list[order(sapply(ff.list, function(x) -length(x)), sapply(ff.list, function(x) x[1]))]  # re-order it nicely
	
	#data.frame(code, numloc, opt.str, SibSize=sapply(ff.list, length), SibSet=sapply(ff.list, paste, collapse=","))	
	ff.list
}

#####################################################
