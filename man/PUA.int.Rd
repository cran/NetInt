% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netint.R
\name{PUA.int}
\alias{PUA.int}
\title{Per-edge Unweighted Average (PUA) network integration}
\usage{
PUA.int(...)
}
\arguments{
\item{...}{a list of numeric matrices. These must be named matrices
representing adjacency matrices of the networks. Matrices may have
different dimensions, but corresponding elements in different matrices
must have the same name.}
}
\value{
the integrated matrix : the matrix resulting from PUA.
}
\description{
It performs the per-edge unweighted average integration between
networks:
\loadmathjax
\mjsdeqn{\bar{w}_{ij} = \frac{1}{|D(i,j)|} \sum_{d \in D(i,j)} w_{ij}^d}
where: \mjsdeqn{D(i,j) = \lbrace  d | v_i \in V^d \wedge v_j \in V^d \rbrace}
}
\examples{
# Create three example networks of different size
set.seed(123);
A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
diag(A1) <- 0;
rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);

A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
diag(A2) <- 0;
rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];

A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
diag(A3) <- 0;
rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));

# Integrate networks using Per-edge Unweighted Average (PUA) method
A_int <- PUA.int(A1, A2, A3);
}
