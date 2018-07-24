\name{VerifyActualParams}
\alias{VerifyActualParams}
\title{
Internal checks to make sure functions are being called with the correctly formatted vector of free parameters.
}
\description{
Silent if all works well, throws an error otherwise.
}

\usage{
VerifyActualParams(x)
}
\arguments{
  \item{x}{vector to check}
   }


\value{
 Ideally nothing, though it throws stop if there is a problem.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
