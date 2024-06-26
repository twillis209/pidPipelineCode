% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h2.R
\name{calcUniQt}
\alias{calcUniQt}
\title{Calculate power to identify heritability of quantitative trait}
\usage{
calcUniQt(n = 1000, hsq = 0.5, alpha = 0.05, var_pi = 2e-05)
}
\arguments{
\item{n}{Sample size}

\item{hsq}{variance explained by all SNPs}

\item{alpha}{significance level}

\item{var_pi}{variance of the off-diagonal elements of the GRM}
}
\value{
list containing se (standard error), ncp (non-centrality parameter) and power
}
\description{
Calculate power to identify heritability of quantitative trait
}
