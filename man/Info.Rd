% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IRF_functions.R
\name{Info}
\alias{Info}
\title{Test information function for 2PL IRF.}
\usage{
Info(parms, theta)
}
\arguments{
\item{parms}{a named list or data.frame with elements \code{parms$alpha} and \code{parms$beta} corresponding to the discrimination and difficulty parameters of the 2PL model, respectively.}

\item{theta}{the latent trait.}

\item{resp}{a matrix or data.frame containing the binary item responses with \code{NA} for missing responses (used omit missing responses).}
}
\value{
\code{length(theta)} vector of test information.
}
\description{
Computes test information under 2PL model.
}
