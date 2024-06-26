% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h2.R
\name{calcBiQt}
\alias{calcBiQt}
\title{Calculate power to identify genetic correlation of quantitative traits}
\usage{
calcBiQt(
  n1 = 1000,
  n2 = 1000,
  hsq1 = 0.5,
  hsq2 = 0.5,
  rg = 0.5,
  rp = 0.5,
  overlap = FALSE,
  alpha = 0.05,
  var_pi = 2e-05
)
}
\arguments{
\item{n1}{sample size of first study}

\item{n2}{sample size of second study}

\item{hsq1}{variance explained by all SNPs in first trait}

\item{hsq2}{variance explained by all SNPs in second trait}

\item{rg}{genetic correlation}

\item{rp}{phenotype correlation}

\item{overlap}{boolean indicating whether there is sample overlap between the two studies}

\item{alpha}{significance level}

\item{var_pi}{variance of the off-diagonal elements of the GRM}
}
\value{
list containing se (standard error), ncp (non-centrality parameter) and power
}
\description{
Calculate power to identify genetic correlation of quantitative traits
}
