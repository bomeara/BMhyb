\name{CheckFlow}
\alias{CheckFlow}
\title{
	CheckFlow
	}

\description{
	Checks to see if the hybridization network is one that is amenable to our algorithm (basically, every lineage may have had no more than one hybridization event leading to it). Returns the vector of recipients that violate this restriction; if the length is zero, then it is a valid network.
	}

\usage{
	CheckFlow(phy,flow)
	}

\arguments{
\item{phy}{your phylogeny}
\item{flow}{a data.frame with gene flow.}
	}



\value{
	A list containing problem.taxa, a vector of recipients that violate the algorithm's restrictions, and all.events, a data.frame of lumped events that lumps into events, rather than sets of taxa.
 	}

\author{
	Dwueng-Chwuan Jhwueng & Brian O'Meara
	}
