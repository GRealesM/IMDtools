## Set of useful functions for file processing and other stuff


##' Convert effect sizes and standard error from linear to log-odds scale
##'
##' Binary traits in case-control studies can be analysed by BOLT-LMM, which uses a linear model.
##' This function provides an approximation from linear to log-odds scale for Beta and SE,
##' as described in BOLT-LMM manual (https://data.broadinstitute.org/alkesgroup/BOLT-LMM/).
##' @title Transform Beta and SE from linear to logg-odds scale
##' @param beta a numeric vector of Beta in linear scale.
##' @param se a numeric vector of SE in linear scale.
##' @param N0 number of controls in the study.
##' @param N1 number of cases in the study.
##' @return a list containing beta and se.beta in log-odds scale.
##' @export
linORscale <- function(beta, se,N0,N1){
	cp  <- N1/(N0+N1)
	BETA  <- beta/(cp * (1-cp))
	SE  <- se/(cp * (1-cp))
	list(BETA,SE)
}

##' Estimate trait standard deviation given vectors of variance of coefficients, MAF and sample size, and correct BETA and SE.
##'
##' For studies of quantitative traits, measurements can come in different scales, leading to
##' differences in standard deviation of the effects.
##' This function estimates the trait standard deviation from MAF and sample size, and by dividing beta and se
##' by this estimate we can make var.beta = 1, thus making different datasets of quantitative traits in different
##' scales comparable.
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##'
##' This is an adaptation of sdY.est function, by Chris Wallace
##'
##' @title Estimate trait variance, internal function
##' @param beta vector of coefficients
##' @param se vector of standard error of coefficients
##' @param maf vector of MAF (same length as beta and se)
##' @param n sample size
##' @return a list containing adjusted Beta and SE by estimated standard deviation of Y
##' @author Guillermo Reales, Chris Wallace
##' @export
sdY.correction <- function(beta, se, maf, n) {
    if(length(beta) != length(maf) | length(se) != length(maf)) stop("Beta, SE and/or MAF are not the same length")
    oneover <- 1/se^2
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    message("Estimated sdY is ", sqrt(cf))
    BETA  <- beta/sqrt(cf)
    SE  <- se/sqrt(cf)
    list(BETA,SE)
}



### g.complement, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)

#' A function to find the complement to a vector of alleles
#'
#' Just like the one in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
#'
#' @param x an allele vector
#'
#' @return a complement vector of alleles
#'
g.complement <- function (x) {
  x <- toupper(x)
  switches <- c(A = "t", T = "a", C = "g", G = "c")
  for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
  toupper(x)
}

#' A function to flip alleles
#'
#' g.rev, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
#'
#' @param x an allele vector
#' @param sep the separator. Default: "/"
#'
#' @return a vector of flipped alleles
#'
g.rev <- function (x, sep = "/") {
  sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
}

#' A function to diagnose differences between two allele vectors
#'
#' g.class, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
#'
#' @param x a first allele vector
#' @param y a second allele vector
#' @param flip_strand If TRUE, it will assume that complementary and reverse complementary alleles are possible. This will also remove ambiguous alleles.
#'
#' @return a diagnostic vector
#'
g.class  <- function(x, y, flip_strand=TRUE){
  diag <- rep(NA, length(x))
  diag[x == g.rev(y)]  <- "rev"
  diag[x == y] <- "nochange"
  if(flip_strand){
    diag[x == g.complement(y)]  <- "comp"
    diag[x == g.rev(g.complement(y))]  <- "revcomp"
    diag[x %in% c("A/T", "T/A", "G/C","C/G")] <- "ambig"
  }
  diag[x != y && (grepl("-", x) | x %in% c("I/D","D/I"))]  <- "indels"
  diag[is.na(diag)]  <- "impossible"
  diag
}

#' A summary statistic allele aligner
#'
#' This function will take a GWAS summary statistic dataset and will allign its alleles to a manifest or reference dataset.
#'
#' GWAS summary statistics should have a minimum number of columns: CHR38, BP38, REF, ALT, BETA, SE, P, and optionally ALT_FREQ.
#' Manifest should have at least CHR38, BP38, REF, ALT.
#'
#' The function will compare both allele placement and check if they're already aligned, if they need to be flipped, find their complement or reverse complementary.
#' Then it will do what's needed for each SNP, including reversing BETA and ALT_FREQ, if necessary.
#'
#'
#' @param ds a data.table containing the GWAS summary statistic dataset
#' @param manifest a data.table containing the manifest, or reference
#'
#' @return a data.table with matching REF and ALT to the reference
#' @export
#'
g.align <- function(ds, manifest){
        # Copy both datasets
		d <- copy(ds)
		man <- copy(manifest)
		
		if(!all(sapply(list(d,man), is.data.table))){
			d <- data.table(d)
			man <- data.table(man)
		}
		# Sanity checks. We'll require these columns initially
		mincold <- c("CHR38", "BP38","REF","ALT", "BETA", "SE", "P")
		if(!all(mincold %in% names(d))) stop("Minimum columns missing from dataset to be aligned, these are: ", paste(mincold, collapse=", "))
		mincolm <- mincold[1:4]
		if(!all(mincolm %in% names(man))) stop("Minimum columns missing from manifest to align dataset to be aligned to, these are: ", paste(mincolm, collapse=", "))

		# Create allele vector for manifest
		man[, alleles := paste(REF,ALT, sep="/")][, pid := paste(CHR38, BP38, sep=":")]
		d[ , REF := toupper(REF)][, ALT:=toupper(ALT)][, alleles:=paste(REF,ALT,sep="/")][, pid := paste(CHR38, BP38, sep=":")]
		M <- merge.data.table(d, man[,.(pid,alleles)], by='pid', suffixes=c(".d",".m"))
        
		# Diagnose alleles
		M[, diag := g.class(alleles.m, alleles.d)]
		if(!all(M$diag == "nochange")){	
	     	    cat("Some SNPs have to be flipped. ", sum(M$diag == "rev"), " to flip, ", sum(M$diag == "comp"), " to find their complement, and ", sum(M$diag == "revcomp"), " to find their reverse complement.\n")
		    M[, alleles.f:= alleles.d]
		    M[diag == "rev", alleles.f:= unlist(g.rev(alleles.d))] # Flip reversed alleles
		    M[diag == "comp", alleles.f:= g.complement(alleles.d)] # Find complement
		    M[diag == "revcomp", alleles.f:= unlist(g.rev(g.complement(alleles.d)))] # Find rev comp
		    
		    # Remove those SNPs that couldn't be successfuly aligned
		    M[, diag2 := g.class(alleles.m, alleles.f)]
		    if(!all(M$diag2 == "nochange")){
	     	    cat("Unfortunately, ", sum(M$diag2 == "ambig"), " SNPs are ambiguous, and ", sum(M$diag2 == "impossible"), " were impossible to align. These will be removed now.\n")
			M <- M[ diag2 == "nochange"]	
		    }
		    
		    # Now update alleles, BETA, and ALT_FREQ
		    M[ diag %in% c("rev", "revcomp"), BETA:= -BETA]
		    if("ALT_FREQ" %in% names(M)){
		    	M[ diag %in% c("rev", "revcomp"), ALT_FREQ:= 1-ALT_FREQ]
		    }
		    M[, c("REF", "ALT") := tstrsplit(alleles.f, "/")] 
	            M[, c("alleles.f",  "diag2") := NULL]
		}

	         M <- unique(M)
             M[, c("alleles.d", "alleles.m", "diag") := NULL] # These columns were included regardless of whether alignment as required
	         if(nrow(M) > nrow(man)) warning("Aligned file has more SNPs than the manifest. Some SNPs might be duplicated.")
	         M
}



