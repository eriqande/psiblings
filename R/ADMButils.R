# R interface to AD Model Builder files
# Author: Hans J. Skaug
# Version: Nov 07.

#' Write to .dat file for admb
#'
#' Writes a named list to a .dat file that admb can read.
#' Basically just write the values of vectors or matrices
#' to the file and puts a comment above each telling which 
#' parameter it is.
#'
#' @param name  path of file to write to
#' @param L list of values to write (ought to be named for clarity)
#' @example: dat_write("epil5sim.dat",list(n=n,p=6,q=q,X=as.matrix(d)))
#' @value Doesn't return anything.  Called for its side effect of writing to a file.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
dat_write <- function(name,L)
{
  n = nchar(name)
  if(substring(name,n-3,n)==".dat")
    file_name = name
  else
    file_name = paste(name,".dat",sep="")

  cat("# \"",name,".dat\" produced by dat_write() from ADMButils; ",date(),"\n", file=file_name,sep="")
  for(i in 1:length(L))
  {
    x = L[[i]]

    if(data.class(x)=="numeric")
      cat("#",names(L)[i],"\n",x,"\n\n",file=file_name,append=T)

    if(data.class(x)=="matrix")
    {
      cat("#",names(L)[i],"\n",file=file_name,append=T)
      write.table(x,,col=F,row=F,quote=F,file=file_name,append=T)
      cat("\n",file=file_name,append=T)
    }

    if(data.class(x)=="list")
    {
      cat("#",names(L)[i],"\n",file=file_name,append=T)
      for(j in 1:length(x))
        if(is.numeric(x[[j]]))
          cat(x[[j]],"\n",file=file_name,append=T)
	else
	  stop("List with non-numeric elements not yet implemented")
      cat("\n",file=file_name,append=T)
    }

  }
}

#' Write to .pin file for admb
#'
#' Writes a named list to a .pin file that admb can read.
#' Probably just write the values of vectors or matrices
#' to the file and puts a comment above each telling which 
#' parameter it is.
#'
#' @param name  path of file to write to
#' @param L list of values to write (ought to be named for clarity)
#' @example: pin_write("kalman_ar1.pin",list(log_sigma=c(0,0),a=0,p=rep(1,10)))
#' @value Doesn't return anything.  Called for its side effect of writing to a file.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
pin_write <- function(name,L)
{
  n = nchar(name)
  if(substring(name,n-3,n)==".pin")
    file_name = name
  else
    file_name = paste(name,".pin",sep="")

  cat("# \"",name,".pin\" produced by pin_write() from ADMButils; ",date(),"\n", file=file_name,sep="")
  for(i in 1:length(L))
  {
    x = L[[i]]

    if(data.class(x)=="numeric")
      cat("#",names(L)[i],"\n",L[[i]],"\n\n",file=file_name,append=T)

    if(data.class(x)=="matrix")
    {
      cat("#",names(L)[i],"\n",file=file_name,append=T)
      write.table(L[[i]],,col=F,row=F,quote=F,file=file_name,append=T)
      cat("\n",file=file_name,append=T)
    }
  }
}


#' Read from a .par file produced by admb
#'
#' Read from a .par file produced by admb or a file of identical structure.
#' This reads the values of the objective function and the named parameters.
#' If some of the parameters are matrices, then their number of columns
#' should be supplied in a named list. (See example)
#' 
#' @param name  path of file to read from
#' @param ncols named list with the number of columns for matrix parameters
#' @example: par_read("sea16.par")
#' @example: par_read("sea16.rep",ncols=list(N=3,ogives=2,age_dist=50)) # N,ogives,age_dist are matrices
#' @value A list with components named as in the .par file, along with three components at the end: $n_par, $loglik, and $gradient which are extracted from the first line of the .par file.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
par_read <- function(name,ncols=list()) # matrices must be specified by name and ncol
{
  n = nchar(name)
  endelse = substring(name,max(1,n-3),n)
  har_endlese = (substring(endelse,1,1) == ".")
  if(har_endlese)
    file_name = name
  else
    file_name = paste(name,".par",sep="")

  tmp = scan(file_name,what="",quiet=T)
  tmp2 = split(tmp,cumsum(tmp=="#"))
  x = tmp2

  if(endelse ==".par")
    x = x[-1]

  for(i in 1:length(x))
  {
    y = x[[i]]
    n = nchar(y[2])
    x[[i]] = as.numeric(y[-(1:2)])
    names(x)[i] = substring(y[2],1,n-1)
  }

  # Convert to matrix for those arguments relevant
  if(length(ncols)>0)
    for(i in 1:length(ncols))
    {
      NN = names(ncols)[i]
      x[[NN]] <- matrix(x[[NN]],ncol=ncols[[i]],byrow=T)
    }

  if(endelse == ".par")
  {
    x$n_par = -as.numeric(tmp2[[1]][6])
    x$loglik = -as.numeric(tmp2[[1]][11])
    x$gradient = -as.numeric(tmp2[[1]][16])
  }

  x
} 


#' Read from a .std file produced by admb
#'
#' Read from a .std file produced by admb or a file of identical structure.
#' 
#' @param name  path of file to read from
#' @example: par_read("sea16.std")
#' @value A list with components est and std for estimates and standard deviations.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
std_read <- function(name)
{
  n = nchar(name)
  if(substring(name,n-3,n)==".std")
    file_name = name
  else
    file_name = paste(name,".std",sep="")

  tmp = read.table(file_name,skip=1)
  est = tmp[,3]
  names(est) = tmp[,2]
  std = tmp[,4]
  names(std) = tmp[,2]

  L1 = list()
  L2 = list()
  
  for(i in unique(names(std)))
  {
    L1[[i]] = est[names(est)==i]
    L2[[i]] = std[names(std)==i]
  }

  list(est=L1,std=L2)
}


#' make a correlation matrix from a covariance matrix
#'
#' If m is a covariance matrix, this returns the corresponding correlation matrix.
#' @param m a covariance matrix
#' @example: par_read("sea16.std")
#' @value A list with components est and std for estimates and standard deviations.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is a similar function in base R or a similar name, but Hans still calls this in places.
cov2corr <- function(m) diag(1/sqrt(diag(m))) %*% m %*% diag(1/sqrt(diag(m)))


#' test each element of x to see if it is in y
#'
#' A function that delivers the same functionality as x %in% y
#'
#' @param x a vector
#' @param y another vector
#' @example: member(1:5, 2:10)
#' @value A logical vector of the same length as x
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is similar functionality in base R, but Hans calls this in places so we include it here.
member <- function(x,y) !is.na(match(x,y))  




#' return logical vector picking out the lower triangle of an n x n square matrix
#'
#' A function that delivers the same functionality as lower.tri(matrix(1,n,n))
#'
#' @param n dimension of a square matrix
#' @param strictly  TRUE means do not include the diagonal and FALSE means do.
#' @example: below(10)
#' @value A logical vector (matrix) with TRUE on the lower diagonal (and on the diagonal with strictly=F)
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is similar functionality in base R, but Hans calls this in places so we include it here.
below <- function(n,strictly=F)
{
  M <- matrix(T,n,n)
  M[rep(1:n,n)<rep(1:n,rep(n,n))] <- F
  if(strictly)
    diag(M) = F
  M
}




#' read Hessian of dimension n from admb .cor file
#'
#' 
#' @param file path to file holding the Hessian
#' @param n number of rows and columns of the Hessian
#' @param cor logical indicating whether to return...
#' @example: readH("boing.cor", 10, T)
#' @value A matrix
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no} 
readH <- function(file,n,cor=F)
{
  N = n*(n+1)/2+4*n
  tmp = scan(file,what="",skip=2,quiet=T)
  if(length(tmp)<N) stop("n is too large")
  tmp = tmp[1:N]

  stdtab = numeric(n)
  H = diag(n)

  for(i in 1:n)
  {
    stdtab[i] = as.numeric(tmp[4])
    tmp = tmp[-(1:4)]
    H[i,1:i] = as.numeric(tmp[1:i])
    tmp = tmp[-(1:i)]
  }

  if(length(tmp)!=0)
  {
    print(length(tmp))
    stop("Det er noe galt")
  }
 
  H = H+t(H)	# Fill in upper diagonal
  diag(H) = 1
  

  if(!cor)
    H = diag(stdtab) %*% H %*% diag(stdtab)

  H
}

