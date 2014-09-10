#' Calls the pseudo likelihood
#' 
#'  Function that fits the pseudo likelihood
#'  @export
pseudo = function(
		A,B,			# Matrices containing SNPs, one individual per line
		position=rep(-1,nrow(A)),# Distance from joining point of stream (-1 indicates missing value)
		stream=rep(1,nrow(A)),	# Stream identifier: 1,2,3
		f,			# Allele frequencies: list of vectors of length 2
					# If missing the calculated empirically from A and B
		eps_error=0,		# Measurement error rate
		relations=c("full-sibs","half-sibs","parent-offspring"), # What should be estimated ?
		fit.regression=all(position>0),# Should regression on "position" be fitted?
		run.model=T,		# If F exe is not called, but files are written/read.
		show.output=F		# If TRUE then output from ADMB is shown on screen
		)
{

  if(!is.matrix(A)) stop("A must be a matrix")
  n = nrow(A)
  S = ncol(A)
  
  # Argument checking 
  if(all(dim(A)!=dim(B))) stop("A and B must be of the same dimension")
  if((!is.vector(position))|(!is.vector(stream))) stop("position and stream must be vectors")
  if((length(position)!=n)|(length(stream)!=n)) stop("position or stream has the wrong length")
  if(!is.numeric(position)) stop("position must be numeric")
  if(!all(member(stream,1:3))) stop("stream can only take the values 1, 2 or 3")

  if(missing(f))
  {
    f = list()
    for(i in 1:S)
    {
      ttt = c(A[,i],B[,i])
      ttt = ttt[ttt>0]	# Remove missing values
      allele_tab = table(ttt)	# Count allele frequencies
      max_allele = max(ttt, na.rm=T)	# Largest allele
      f[[i]] = numeric(max_allele)
      names(f[[i]]) = 1:(max_allele)
      f[[i]][names(allele_tab)] = allele_tab/length(ttt)
    }
  }
  
  # Find largest allele at all loci
  max_allele = unlist(lapply(f,length))

  # Write to file-----------------------

  # Data and starting values  
  dat_write("pseudo2",L=list(n=n,S=S,A=A,B=B,
  	position=position,stream=stream,max_allele=max_allele,f=f,
  	eps_error=eps_error))
  #
  # Control statements  	
  cat("# Fit regression coeffs (1) or not (0)","\n",as.numeric(fit.regression),"\n",file="phase_pseudo.dat") 
  cat("# Fit parent offspring (1) or not (0)","\n",as.numeric(member("parent-offspring",relations)),
  		"\n",file="phase_pseudo.dat",append=T) 
  cat("# Print posterior probabilties (1) or not (0)","\n",1,"\n",file="phase_pseudo.dat",append=T)
  
  # Fit model and read back results ----------------------------
  if(run.model)
  {
    system(paste(Pseudo2BinaryPath(), " -gbs 3000000000 -est"),ignore.stdout=!show.output)		# Run executable
    fit = par_read("pseudo2.par")
    fit$n = n
    fit$A = A
    fit$B = B
    fit$f = f
    fit$relations = relations
    fit$eps_error = eps_error
    fit$posterior = read.table("pseudo2.rep")	# Read posterior probabilities
	    colnames(fit$posterior) = c("i","j","Unrelated","Full sibs","Half sibs","Parent-offspring")
    tmp = fit$posterior[,"Full sibs"]
    fit$posterior$logit_P_fs = log(tmp/(1+1e-10-tmp))
    fit$relations = relations
    return(fit)
  }
  else
    return(invisible())
}

# Determine set of full-sibling groups
full_siblings = function(fit,c.tuning=1,nsim=1)
{
  # Setting c.tuning<1 make	s the method less agreesive
  # If c.tuning is a vector, the element minimizing a simulation based error rate is chosen (see below)
  # nsim = number of times simulations are repeated to evaluate error rate

  if(length(unique(c.tuning))>length(c.tuning))
    stop("Elements in \"c.tuning\" must be unique")
#  if(any((c.tuning<=0)|(c.tuning>1)))
#    stop("\"c.tuning\" can only take values in (0,1)")

  # Packages
  require(igraph)

  n = fit$n					# Sample size
  N = n*(n-1)/2					# Number of pairwise comparisons
  p_fs = exp(fit$a)/(1+exp(fit$a))  		# Estimated fraction of full sib pairs 
  Mpost = fit$posterior				# Matrix of posterior probabilities

  
  # Må låse sannsynligheter for at dette skal virke
  # Try to testimate the true number of pairs in the dataset
  # Summing posterior probabilities, and subtracting the corresponding values under the null hypothesis 
  # of no related individuals
  if(F)
  {
      # All unrelated
      Apert = rowwise_perturb(fit$A)
      Bpert = rowwise_perturb(fit$B)
      fit0 = pseudo(A=Apert,B=Bpert,relations=fit$relations,eps_error=fit$eps_error,show.output=F)
      n_hat = (sum(fit$posterior[,4]-fit0$posterior[,4]))
      print(sum(fit0$posterior[,4]))
      print(n_hat)
  }

  
  # Selects the top p_fs-fraction of list sorted according to P_posterior(full sib) 
  #kk=Mpost[Mpost[,4] > quantile(Mpost[,4],1.0-c.tuning[1]*p_fs),]
  #	n1 = nrow(kk)

  n1=sum(Mpost[,4])
  kk=Mpost[Mpost[,4] > quantile(Mpost[,4],1.0-c.tuning[1]*n1/N),]
  n1 = nrow(kk)


  if(n1==0)
    output = NULL				# No full sibling groups found
  else
  {
    navn = rownames(fit$A)
    g = graph.data.frame(data.frame(navn[kk[,1]],navn[kk[,2]]),directed=F)
    #E(g)$weight = pmax(kk$logit_P_fs,1e-10)
   
    #f1=leading.eigenvector.community(g)			# Clustering algorithm
    #f1=edge.betweenness.community(g)			# Clustering algorithm
    #f1=fastgreedy.community(g,modularity=T,merges=T)			# Clustering algorithm

     f1 = clusters(g)
    output = split(V(g)$name,f1$mem)
  }

  
  # Find optimal "c.tuning" by simulation 
  if(length(c.tuning)>1)
  {        
    # Determine the set of full and half sibs that should not be perturbed
    full_siblings = unique(c(kk[,1],kk[,2]))
    p_hs = exp(fit$a_hs)/(1+exp(fit$a_hs))  		# Estimated fraction of full sib pairs 
    kk_hs = Mpost[Mpost[,5] > quantile(Mpost[,5],1.0-c.tuning[1]*p_hs),]
    half_siblings = unique(c(kk_hs[,1],kk_hs[,2]))
    siblings = unique(c(full_siblings,half_siblings))
    ind_pert = !member(1:n,siblings)

    # Vector in which "distances are stored"
    dist_tab = numeric(length(c.tuning))    

    for(K in 1:nsim)
    {
      # Perturbed genotypes
      Apert = fit$A
      Bpert = fit$B
    
      # Permutation of alleles for all individuals who is neither half or full sibling
      Apert[ind_pert,] = rowwise_perturb(Apert[ind_pert,])
      Bpert[ind_pert,] = rowwise_perturb(Bpert[ind_pert,])
      if(K==1)
      {
        cat(paste("\n\n-----\nfull_siblings: proportion of individuals perturbed: ",sum(ind_pert)/n,"\n"))
        cat(paste("full_siblings: number of full siblings kept in simulations: ",length(full_siblings),"\n"))
        cat(paste("full_siblings: number of half siblings kept in simulations: ",length(half_siblings),"\n"))
      }

      # Permute also half siblings, but genotype-wise to "keep" the half-sib relationship
      for(k in 1:nrow(kk_hs))
      {
        i = kk_hs[k,"i"]	
        j = kk_hs[k,"j"]	
        tmpind = sample(c(T,F),size=ncol(Apert),replace=T)	# Indicator (by locus) for whehter or not to permute "i" and "j" 
        tmpAi = Apert[i,]						# Take backup
        tmpBi = Bpert[i,]
        Apert[i,tmpind] = Apert[j,tmpind]
        Bpert[i,tmpind] = Bpert[j,tmpind]
        Apert[j,tmpind] = tmpAi[tmpind]
        Bpert[j,tmpind] = tmpBi[tmpind]
      }
    
      outputlist = list()
      # Run through all possible c.tuning values and calculate error
      for(i in 1:length(c.tuning))
      {
        #browser()
        fit_i = pseudo(A=Apert,B=Bpert,relations=fit$relations,eps_error=fit$eps_error,show.output=F)
        outputlist[[i]] = full_siblings(fit=fit_i,c.tuning=c.tuning[i])
        dist_tab[i] = dist_tab[i] + naive_metric(output,outputlist[i])
      }
    
    } # Simulation loop
    
    # Find value of i gives the smalles smallest distance. IN case there is more than 1 such i-value, take largest c.tune
    i_min = which(dist_tab==min(dist_tab))
    i_min = i_min[which(c.tuning[i_min]==max(c.tuning[i_min]))]
    
    # Rerun with original data but optimal c.tuning
    output = full_siblings(fit=fit,c.tuning=c.tuning[i_min])
    cat(paste("full_siblings: tuning parameter chosen: ",c.tuning[i_min]," (distance=",dist_tab[i_min],")\n\n"))

  } #if(length(c.tuning)>1)

  return(output)
}

# perturb allele matrices rowise
rowwise_perturb = function(M)
{
  for(i in 1:ncol(M))
    M[,i] = sample(M[,i])
  return(M)
}

# ---------- Metric used for distance between true and inferred sibroups
 
# Simply the number of individuals not common to the two lists (regardsless of sibling structure)
# Gives more weight to false sibgroups
naive_metric = function(L_true,L_estimated)
{
  ind1 = unlist(L_true)
  ind2 = unlist(L_estimated)
  return(sum(!member(ind1,ind2))+sum(!member(ind2,ind1)))
}
  


