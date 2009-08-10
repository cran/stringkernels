
# wordkernel
setClass("stringkernelEx", representation(name="character",preprocessor="function"), contains="stringkernel")


worddot <- function(type = c("spectrum", "boundrange", "constant","exponential"),
    length = 4, lambda = 1.1, normalized=TRUE, tokenizer=openNLP::tokenize)
{
  type = match.arg(type)
  ktypes = c(constant=1, exponential=2, spectrum=3, boundrange=4)

  if (type=="exponential") param = lambda
  else param = length
  
  rval <- function(xtxt,ytxts=NULL, normalized=FALSE, preprocess=TRUE) {
    if(preprocess) {
      xtxt = tokenizer(xtxt)
      if (!is.null(ytxts)) ytxts = lapply(ytxts, tokenizer)
    }
    
    if (normalized==TRUE) {
      if (is.null(ytxts)) return(1)
      else if (is.list(ytxts) && length(ytxts)!=1) 
        stop("Normalization not supported for lists of y.\n Use kernelMatrix instead")
      else { 
        return(rval(xtxt,ytxts,normalized=FALSE,preprocess=FALSE)/
              sqrt(rval(xtxt,NULL,normalized=FALSE,preprocess=FALSE)*
              rval(unlist(ytxts),NULL,normalized=FALSE,preprocess=FALSE)))
      }
    }

    

    alphabet = unique(xtxt)
    ## index
    idx = structure(32+1:length(alphabet), names=alphabet)
    x = idx[xtxt]
    if (is.null(ytxts)) {
      y = list(x)
    }
    else {
      if (!is.list(ytxts)) ytxts = list(ytxts)
      y = lapply(ytxts, function(u) {v = idx[u]; v[is.na(v)] = 32; v})
    }
    yext = lapply(y, c, 10)
    yext = lapply(yext, as.integer)
    leny = sapply(yext, length)
    

    return(.Call("wordkernel", as.integer(c(x,10)), yext,
              as.integer(length(yext)), as.integer(length(x)+1), as.integer(leny),
              as.integer(ktypes[type]), as.double(param)))
  }
  

  return(new("stringkernelEx", .Data = rval, name="Worddot", kpar = list(length = length,
        lambda = lambda, type = type, normalized = TRUE), preprocessor=tokenizer))
}


stringdotEx <- function(type = c("spectrum", "boundrange", "constant","exponential"),
    length = 4, lambda = 1.1, normalized=TRUE) 
{
  type = match.arg(type)
  ktypes = c(constant=1, exponential=2, spectrum=3, boundrange=4)

  if (type=="exponential") param = lambda
  else param = length
  
  rval <- function(xtxt,ytxts=NULL, normalized=FALSE, preprocess=FALSE) {
    
    if (normalized==TRUE) {
      if (is.null(ytxts)) return(1)
      else if (is.list(ytxts) && length(ytxts)!=1) 
        stop("Normalization not supported for lists of y.\n Use kernelMatrix instead")
      else { 
        return(rval(xtxt,ytxts,normalized=FALSE,preprocess=FALSE)/
              sqrt(rval(xtxt,NULL,normalized=FALSE,preprocess=FALSE)*
              rval(unlist(ytxts),NULL,normalized=FALSE,preprocess=FALSE)))
      }
    }

    # concatenate texts, replace \n with \r and add \n as sentinel
    prep = function(x) paste(paste(chartr("\n", "\r", x), collapse="\r"), "\n", sep="")
    x = prep(xtxt)
    if (is.null(ytxts)) {
      y = list(x)
    } else { 
      y = sapply(ytxts, prep)
    }


    return(.Call("stringtv", as.character(x), as.character(y),
              as.integer(length(y)), as.integer(nchar(x)), as.integer(nchar(y)),
              as.integer(ktypes[type]), as.double(param)))
  }
  

  return(new("stringkernelEx", .Data = rval, name="stringdot", kpar = list(length = length,
        lambda = lambda, type = type, normalized = TRUE), preprocessor=identity))
}


gapweightkernel <- function(length=2, lambda=0.75, 
      normalized=TRUE, tokenizer=openNLP::tokenize, use_characters=FALSE) 
{ 
  if (use_characters) tokenizer = function(x) strsplit(x, "")[[1]]
  obj = new("stringkernelEx", .Data = stop, name="gapweighted", kpar = list(length = length,
        lambda = lambda, normalized = normalized), preprocessor=tokenizer)
  
  rval <- function(xtxt,ytxts=NULL, normalized=obj@kpar$normalized, preprocess=TRUE) {
    if(preprocess) {
      xtxt = obj@preprocessor(xtxt)
      if (!is.null(ytxts)) ytxts = lapply(ytxts, obj@preprocessor)
    }


    if (normalized==TRUE) {
      if (is.null(ytxts)) return(1)
      else if (is.list(ytxts) && length(ytxts)!=1) 
        stop("Normalization not supported for lists of y.\n Use kernelMatrix instead")
      else { 
        return(rval(xtxt,ytxts,normalized=FALSE, preprocess=FALSE)/
              sqrt(rval(xtxt,NULL,normalized=FALSE, preprocess=FALSE)*
              rval(unlist(ytxts),NULL,normalized=FALSE, preprocess=FALSE)))
      }
    }
    
    
    alphabet = unique(xtxt)
    ## index
    idx = structure(1:length(alphabet), names=alphabet)
    x = idx[xtxt]
    if (is.null(ytxts)) {
      y = list(x)
    }
    else {
      if (!is.list(ytxts)) ytxts = list(ytxts)
      y = lapply(ytxts, function(u) {v = idx[u]; v[is.na(v)] = 0; v})
    }
    y = lapply(y, as.integer)
    leny = sapply(y, length)
    
    ret = c();
    for (i in 1:length(y)) {
      ret[i] = .Call("gapkernel_range", as.integer(x), as.integer(y[[i]]), as.integer(length(x)), 
        as.integer(length(y[[i]])), as.integer(max(x)+1), as.integer(obj@kpar$length), 
         as.integer(obj@kpar$length), as.double(obj@kpar$lambda))
    }
    return(ret)
  }
  obj@.Data = rval
  return(obj)
  
}

setMethod("kernelMatrix", 
    signature="stringkernelEx", 
    function(kernel, x, y=NULL) 
{
  if (!is.list(x) || (!is.null(y) && !is.list(y))) 
    stop("kernelMatrix for wordkernel accepts only lists")
  
  x = lapply(x, kernel@preprocessor)
  lx = length(x)
  
  # if y is null, create a symmetric Matrix K_ij = {k(x_i,x_j)}
  if (is.null(y)) {
    res = matrix(0, lx, lx)
    for (i in 1:lx) {
      res[i, i:lx] = kernel(x[[i]], x[i:lx], normalized=FALSE, preprocess=FALSE)
    }
    d = diag(res)
    if (kpar(kernel)$normalized==TRUE) {
      for (i in 1:lx) {
        res[i,i:lx] = res[i,i:lx] / sqrt(d[i] * d[i:lx])
      }
      res = res + t(res)
      diag(res) = 1
    }
    else {
      res = res + t(res)
      diag(res) = d
    }
  }
  else {
    y = lapply(y, kernel@preprocessor)
    ly = length(y)
    res = matrix(0,lx,ly)
    for (i in 1:lx) {
      res[i,] = kernel(x[[i]], y, normalized=FALSE, preprocess=FALSE)
    }
    if (kpar(kernel)$normalized==TRUE) {
      normx = sapply(x, kernel, normalized=FALSE, preprocess=FALSE)
      normy = sapply(y, kernel, normalized=FALSE, preprocess=FALSE)
      res = res / sqrt(normx%o%normy)
    }
  }
  return(as.kernelMatrix(res))
}
)

setMethod("kernelMult", "stringkernelEx",
          function(kernel, x, y=NULL, z, blocksize = 256)
{
  return(kernelMatrix(kernel,x,y)%*%z)
}
)
    
setMethod("kernelPol", "stringkernelEx",
    function(kernel, x, y=NULL ,z ,k=NULL) stop("kernelPol is not supported yet for stringkernelEx"))

setMethod("show", "stringkernelEx", 
    function(object) {cat("String-based kernel", object@name, "\n")})


setClass("stringkernelPrecomputed", representation(items="list", cache="matrix", 
  use_dummy="logical"), contains="stringkernelEx")

precomputedkernel <- function(kernel, items, kernelmatrix=NULL, use_kernel=FALSE, use_dummy=TRUE) {
  
  if (is.character(kernel)) {
    if (is.null(kernelmatrix) || use_kernel) 
      stop("You have to supply a kernelmatrix and use_kernel must be FALSE")
    kernel = new("stringkernelEx", .Data = stop, name=kernel, 
      kpar = list(), preprocessor=identity)
  }
  
  
  if (!is(kernel, "stringkernelEx")) stop("kernel needs to be of class stringkernelEx")

  if (is.null(kernelmatrix)) kernelmatrix = kernelMatrix(kernel, items)
    
  obj = new("stringkernelPrecomputed", .Data = stop, name=paste(kernel@name, "(precomputed)"), 
  kpar = kernel@kpar, preprocessor=kernel@preprocessor, items=items, cache=kernelmatrix, use_dummy=use_dummy)
  
  if (use_kernel) {
    kfunc = kernel@.Data
  } else {
    kfunc <- function(...) stop("use of kernel function not supported")
  }
  
  listmatch <- function(x,l) {
    for(i in 1:length(l))
    if (identical(x,l[[i]])) return(i)
    return(NA)
  }
  
  func <- function(xtxt,ytxts=NULL, normalized=kernel@normalized, preprocess=TRUE) {
    xpos = listmatch(xtxt,obj@items)
    ypos = sapply(ytxts, listmatch, l=obj@items)
    if (is.na(xpos)) return(kfunc(xtxt,ytxts,normalized,preprocess))
    if (any(is.na(ypos))) return(kfunc(xtxt,ytxts,normalized,preprocess)) # TODO improve
    
    return(obj@cache[xpos,ypos])
  } 
  obj@.Data = func
  return(obj)
}

setMethod("kernelMatrix", "stringkernelPrecomputed",
  function(kernel, x, y=NULL) {
  
    listmatch <- function(x,l) {
      for(i in 1:length(l))
      if (identical(x,l[[i]])) return(i)
      return(NA)
    }
    if (!is.list(x) || (!is.null(y) && !is.list(y))) 
      stop("kernelMatrix for stringkernel accepts only lists")
    
    xpos = if (kernel@use_dummy) as.integer(x) else sapply(x, listmatch, l=kernel@items)
    if (is.null(y)) {
      if(any(is.na(xpos))) stop("Item not found in cache!", x[is.na(xpos)])
      return(as.kernelMatrix(kernel@cache[xpos,xpos]))
    } else {
      ypos = if (kernel@use_dummy) as.integer(y) else sapply(y, listmatch, l=kernel@items)
      if(any(is.na(ypos))) stop("Item not found in cache!",  y[is.na(ypos)])
      return(as.kernelMatrix(kernel@cache[xpos,ypos]))
    }
  }
)


setClass("multigapweight", representation(.Data="array", items="list", 
  minlength="vector", maxlength="vector", kpar="list", tokenizer="function"))


multigapweightkernel <- function(items, maxlength, kernelarray=NULL, lambda=0.75, 
      normalized=TRUE, tokenizer=openNLP::tokenize, minlength=1) 
{
  if (is.null(kernelarray)) {
    len = length(items)
    span = maxlength-minlength+1
    kernelarray = array(0, dim=c(len,len,span))
    
    itemsp = lapply(items, tokenizer)
    for (i in 1:len) {
      alphabet = unique(itemsp[[i]])
      ## index
      idx = structure(1:length(alphabet), names=alphabet)
      x = idx[itemsp[[i]]]
      y = lapply(itemsp[i:len], function(u) {v = idx[u]; v[is.na(v)] = 0; v})
      
      for (j in 1:(len-i+1)) {
        kval = .Call("gapkernel_range", as.integer(x), as.integer(y[[j]]), as.integer(length(x)), 
          as.integer(length(y[[j]])), as.integer(max(x)+1), as.integer(maxlength), 
          as.integer(minlength), as.double(lambda))
        kernelarray[i,i+j-1,1:span] = kval
      }
    }
  }

  return(new("multigapweight", .Data=kernelarray, items=items, minlength=minlength, 
      maxlength=maxlength, kpar=list(lambda=lambda, normalized=normalized), 
      tokenizer=tokenizer))
}

setGeneric("getkernel", function(mgw, length, ...) standardGeneric("getkernel"))
setMethod("getkernel", "multigapweight", 
  function(mgw, length, use_dummy=FALSE)
{
    i = length - mgw@minlength + 1
    km =  mgw@.Data[,,i]
    
    
    dkm = diag(km)
    
    
    if (mgw@kpar$normalized) {
      ds = sqrt(diag(km))
      km = km + t(km)
      km = km / (ds%o%ds)
      diag(km) = 1
    } else {
      dkm = diag(km)
      km = km + t(km)
      diag(lm) = dkm
    }
    
    return(precomputedkernel(gapweightkernel(length, mgw@kpar$lambda, mgw@kpar$normalized, mgw@tokenizer),
      mgw@items, as.kernelMatrix(km), use_dummy=use_dummy))
}
)    

precomputeddummy <- function(items) lapply(1:length(items), as.character)
