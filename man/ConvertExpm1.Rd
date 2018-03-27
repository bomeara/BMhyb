\name{ConvertExpm1}
\alias{ConvertExpm1}
\title{
	Exponentials
	}

\description{
	'ConvertExpm1' computes exp(x) - 1 accurately also for |x| << 1.
	}

\usage{
	ConvertExpm1(x)
	}

\arguments{
  \item{x}{
 x: a numeric or complex vector
 }
}

\details{
  see function expm1 in for more detail descrition.
	}

\value{
 A vector of the same length as 'x' containing the transformed value.
 }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. 1988. The New S Language.  Wadsworth & Brooks/Cole.  (for ‘log’, ‘log10’ and ‘exp’.)
	}
