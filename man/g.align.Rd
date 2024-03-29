% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gr-wrenches.R
\name{g.align}
\alias{g.align}
\title{A summary statistic allele aligner}
\usage{
g.align(ds, manifest)
}
\arguments{
\item{ds}{a data.table containing the GWAS summary statistic dataset}

\item{manifest}{a data.table containing the manifest, or reference}
}
\value{
a data.table with matching REF and ALT to the reference
}
\description{
This function will take a GWAS summary statistic dataset and will allign its alleles to a manifest or reference dataset.
}
\details{
GWAS summary statistics should have a minimum number of columns: CHR38, BP38, REF, ALT, BETA, SE, P, and optionally ALT_FREQ.
Manifest should have at least CHR38, BP38, REF, ALT.

The function will compare both allele placement and check if they're already aligned, if they need to be flipped, find their complement or reverse complementary.
Then it will do what's needed for each SNP, including reversing BETA and ALT_FREQ, if necessary.
}
