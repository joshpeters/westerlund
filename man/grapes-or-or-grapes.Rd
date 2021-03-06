% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{\%||\%}
\alias{\%||\%}
\title{Set a default value if an object is null}
\usage{
lhs \%||\% rhs
}
\arguments{
\item{lhs}{An object to set if it's null}

\item{rhs}{The value to provide if x is null}
}
\value{
rhs if lhs is null, else lhs
}
\description{
Set a default value if an object is null
}
\references{
https://adv-r.hadley.nz/functions.html#missing-arguments
}
\author{
Hadley Wickham
}
