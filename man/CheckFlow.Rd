\name{CheckFlow}
\alias{CheckFlow}
\title{
	CheckFlow
	}

\description{
	Checks to see if the hybridization network is one that is amenable to our algorithm (basically, every lineage may have had no more than one hybridization event leading to it). Returns the vector of recipients that violate this restriction; if the length is zero, then it is a valid network.
	}

\usage{
	CheckFlow(flow)
	}

\arguments{
\item{flow}{a data.frame with gene flow.}
	}



\value{
	Vector of recipients that violate the algorithm's restrictions
 	}

\author{
	Dwueng-Chwuan Jhwueng & Brian O'Meara
	}
