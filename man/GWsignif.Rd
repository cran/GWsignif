\name{GWsignif}
\alias{GWsignif}
\title{
Estimating genome-wide significance for whole genome sequencing studies
}
\description{
Extrapolate whole genomic region significance threshold from a small region significance threshold. 
}
\usage{
GWsignif(pvalues, ngroup = 5, ntest.genome, alpha = 0.05, plot.figure = TRUE)
}
\arguments{
  \item{pvalues}{ a matrix or data frame of permutation pvalues in sub-regions, or a vector of files storing the sub-regions pvalues, no header in files.}
  \item{ngroup}{ the number of groups. The sub-regions of same size are put in a same group.}
  \item{ntest.genome}{ the number of tests in a large genome-wide region of interest.}
  \item{alpha}{ family-wise error rate (FWER), default 0.05.}
  \item{plot.figure}{ plot a figure if TRUE.}
}
\details{
Estimate genome-wide significance using extrapolation approach: estimate significance thresholds in small regions and then extrapolate a large region significance threshold.
Procedure of the extrapolation approach:
1) Define a small region of size m along Chromosomes and divide the samll region into sub-regions of size m/2^(k-1), k=1, 2, 3, ..., K. It is better to choose the small region appropriately large, but small compared to whole genome.  
2) Perform permutation or simulation tests for each sub-region: 
   -  repeat each test in each region using normally-distributed phenotypes (under the null) or permutations of the given phenotypes;
   -  calculate region significance from the repeated tests in each region by alpha-quantile (e.g., 5\%) of the minimum p-values in the region.
3) Model the region significance as a function of the region size using simple linear regression.
4) Predict (extrapolate) the significance threshold in a large whole genomic region using the model.
}
\value{
  \item{qminp }{ a list of alpha-quantile of minimum pvalues within a sub-region in different groups, that is, empirical significance threshold in each sub-region.}
  \item{mlogq }{ a matrix consisting of ntest (number of tests in a group), bonf (-log10 of Bonferroni correction, i.e., -log10(alpha/ntest)), mean and standard error of -log10(qminp). The row represents a different group.}
  \item{alpha }{ a desired family-wise error rate (FWER)}
  \item{ntest.region }{ the number of tests in the small region defined before.}
  \item{ntest.genome }{ the number of tests in the large genome-wide region of interest}
  \item{GWsignif.threshold }{ estimated significance threshold in the large genome-wide region}
}
\references{
ChangJiang Xu, Ioanna Tachmazidou, Klaudia Walter, Antonio Ciampi, Eleftheria Zeggini, Celia M.T. Greenwood (2014) Estimating genome-wide significance for whole genome sequencing studies. Genetic Epidemiology
}
\author{ChangJiang Xu}
\examples{
##Example
##pvaluesDir <- "pvalues/"
##pvalues <- paste0(pvaluesDir, pvalues_file_names)
##gw <- GWsignif(pvalues, alpha = 0.05)

set.seed(630657)
pvalues <- matrix(runif(1000*5000), 1000)
gw <- GWsignif(pvalues, alpha = 0.05)
gw <- GWsignif(pvalues, ngroup = 5, ntest.genome = 50000, alpha = 0.05)
}
