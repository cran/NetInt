# Network Integration library
# November 2012
# January 2013
# October 2013
# September 2022 - Created package for CRAN: NetInt 1.0.0

###############################################################################
################## Unweighted Network Integration #############################
###############################################################################

#' Unweighted Average (UA) network integration
#'
#' @description It performs the unweighted average integration between networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij} = \frac{1}{n} \sum_{d = 1}^n w_{ij}^d}
#'
#' @param align logical. If TRUE (def.) the matrices are aligned using
#'      align.networks, otherwise they are directly summed without any
#'      previous alignment.
#' @param ... a list of numeric matrices. These must be named matrices
#'       representing adjacency matrices of the networks. Matrices may have
#'       different dimensions, but corresponding elements in different matrices
#'       must have the same name.
#'
#' @return the integrated matrix : the matrix resulting from UA.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Integrate networks using Unweighted Average (UA) method
#' A_int <- UA.int(align=TRUE, A1, A2, A3);
UA.int <- function(align=TRUE, ...) {
    if (align)
        nets <- align.networks(fill=0, ...)
    else
        nets <- list(...)
    res <- Reduce("+", nets)
    return(res/length(nets))
}



#' "Memory Saving" Unweighted Average (UA) network integration
#'
#' @description It performs the unweighted average integration between networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij} = \frac{1}{n} \sum_{d = 1}^n w_{ij}^d}
#' The matrices are read from files and loaded one at time in memory.
#'
#' @param nets.files a list with the names of the .rda files storing the
#'      matrices representing the weighted adjacency matrices of the nets.
#' @param example.names a list with the names of the examples stored in each
#'      net. If NULL (def.) it is assumed that the matrices have exactly the
#'      same examples in the same order, otherwise the matrices are arranged to
#'      be correctly aligned.
#'
#' @return the integrated matrix : the matrix resulting from UA.
#' @export
MS.UA.int <- function(nets.files, example.names=NULL) {
    n <- length(nets.files)
    if (is.null(example.names)) {
        data.name <- load(nets.files[[1]])
        M <- eval(parse(text=data.name))  # M represents the adjacency matrix
        if (n>1)
            for(i in 2:n) {
                data.name <- load(nets.files[[i]])
                N <- eval(parse(text=data.name))
                M <- M + N
                rm(N)
                gc()
            }
    } else {
        ex <- unique(unlist(example.names))
        m <- length(ex)
        M <- matrix(numeric(m*m),nrow=m)
        rownames(M) <- colnames(M) <- ex
        for(i in 1:n) {
            data.name <- load(nets.files[[i]])
            N <- eval(parse(text=data.name))
            exN <- example.names[[i]]
            M[exN,exN] <- M[exN,exN] + N
            rm(N)
            gc()
        }
    }
    return(M/n)
}



#' Per-edge Unweighted Average (PUA) network integration
#'
#' @description It performs the per-edge unweighted average integration between
#' networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij} = \frac{1}{|D(i,j)|} \sum_{d \in D(i,j)} w_{ij}^d}
#' where: \mjsdeqn{D(i,j) = \lbrace  d | v_i \in V^d \wedge v_j \in V^d \rbrace}
#'
#' @param ... a list of numeric matrices. These must be named matrices
#'       representing adjacency matrices of the networks. Matrices may have
#'       different dimensions, but corresponding elements in different matrices
#'       must have the same name.
#'
#' @return the integrated matrix : the matrix resulting from PUA.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Integrate networks using Per-edge Unweighted Average (PUA) method
#' A_int <- PUA.int(A1, A2, A3);
PUA.int <- function(...) {
    dummy.value <- 10^10
    n <- length(list(...))

    # counting missed values for each network
    nets<-align.networks(fill=dummy.value, ...)

    for (i in 1:n) {
        m <- nets[[i]]
        m[m!=dummy.value] <- 0
        m[m==dummy.value] <- 1
        nets[[i]] <- m
    }
    nets.missed <- Reduce("+", nets)

    # summing values for each network
    nets <- align.networks(fill=0, ...)
    nets.sum <- Reduce("+", nets)

    # normalizing for the number of edges effectively present
    nets.missed[nets.missed==n] <- 1 #this is to avoid division by zero. Note that in this case nets.sum==0
    return(nets.sum/(n-nets.missed))
}



#' Maximum (MAX) network integration
#'
#' @description It performs the Max integration between networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij} = \max_{d} w_{ij}^d}
#'
#' @param ... a list of numeric matrices. These must be named matrices
#'       representing adjacency matrices of the networks. Matrices may have
#'       different dimensions, but corresponding elements in different matrices
#'       must have the same name.
#'
#' @return the integrated matrix : the matrix resulting from MAX.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Integrate networks using Maximum (MAX) method
#' A_int <- MAX.int(A1, A2, A3);
MAX.int <- function(...) {
    nets<-align.networks(fill=0, ...)
    res <- Reduce("pmax", nets)
    return(res)
}



#' Minimum (MIN) network integration
#'
#' @description It performs the Min integration between networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij} = \min_{d} w_{ij}^d}
#' Note that this function consider the minimum between existing edges, that is
#' if an edge (i,j) is not present in a network, since one of the nodes
#' i or j is not present in the network, then the edge is not considered in the
#' computation. If the edge (i,j) is not present in any of the available
#' networks, that its value is 0. If drastic=TRUE the minimum is zero if at least
#' one edge is not present in a network.
#'
#' @param drastic if TRUE the minimum is zero if at least one edge is not present
#'          in a network (def: FALSE).
#' @param ... a list of numeric matrices. These must be named matrices representing
#'       adjacency matrices of the networks. Matrices may have different dimensions,
#'       but corresponding elements in different matrices must have the same name.
#'
#' @return the integrated matrix : the matrix resulting from MIN.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Integrate networks using Minimum (MIN) method
#' A_int.noDrastic <- MIN.int(drastic = FALSE, A1, A2, A3);
#'
#' # Integrate networks using Minimum (MIN) method (drastic integration)
#' A_int.drastic <- MIN.int(drastic = TRUE, A1, A2, A3);
MIN.int <- function(drastic=FALSE, ...) {
    dummy.max <- 10000
    if (drastic)
        filled <- 0
    else
        filled <- dummy.max
    nets<-align.networks(fill=0, ...)
    n.nets <- length(nets)
    if (!drastic)
        for (i in 1:n.nets)
            nets[[i]][nets[[i]] == 0] <- dummy.max
    res <- Reduce("pmin", nets)
    if (!drastic)
        res[res==dummy.max] <- 0
    return(res)
}



#' ATLEASTK network integration
#'
#' @description It performs the ATLEAST integration between networks: only edges
#' present in at least k networks are preserved, the others are eliminated.
#' The resulting network is a binary network: the edge is 1 if preserved, otherwise 0.
#' An edge is considered "present" if its value is larger than 0.
#'
#' @param k the minimum number of the networks in which each edge must be present to
#'     be preserved (k=1).
#' @param ... a list of numeric matrices. These must be named matrices representing
#'       adjacency matrices of the networks. Matrices may have different dimensions,
#'       but corresponding elements in different matrices must have the same name.
#'
#' @return the integrated matrix : the matrix resulting from ATLEASTK.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Integrate networks using ATLEASTK method
#' A_int <- ATLEASTK.int(k=2, A1, A2, A3);
ATLEASTK.int <- function(k=1, ...) {
    n <- length(list(...))

    # counting missed values for each network
    nets<-align.networks(fill=0, ...)

    for (i in 1:n) {
        m <- nets[[i]]
        m[m>0] <- 1
        m[m<=0] <- 0
        nets[[i]] <- m
    }
    nets.counts <- Reduce("+", nets)
    nets.counts[nets.counts<k] <- 0
    nets.counts[nets.counts>=k] <- 1

    return(nets.counts)
}



###############################################################################
################## Weighted Network Integration ###############################
###############################################################################

#' Weighted Average Per-class (WAP) network integration
#'
#' @description It performs the WAP integration between networks:
#' \loadmathjax
#' \mjsdeqn{\bar{w}_{ij}(k) = \sum_{d = 1}^n \alpha^d(k) w_{ij}^d}
#' where
#' \mjsdeqn{\alpha^d(k) = \frac{1}{\sum_{j=1}^n M^j(k)} M^d(k)}
#' and \mjseqn{M^d(k)} is a suitable accuracy metrics for class k on network d.
#' The metrics could be, e.g. the AUC or the precision at a given recall.
#' Note that this function puts more weight (alpha parameter) for networks with
#' associated larger M.
#'
#' @param m a numeric vector with the values of the metric used to compute the alpha
#'     coefficients. It could be e.g. AUC values.
#' @param align logic. If TRUE the numeric matrices passed as arguments are aligned
#'         according to the function align.networks (def: FALSE).
#' @param logint logic. If TRUE m is log transformed: -log(1-m), otherwise a linear
#'          integration is performed (def: FALSE).
#' @param ... a list of numeric matrices. These must be named matrices representing
#'       adjacency matrices of the networks. Matrices may have different dimensions,
#'       but corresponding elements in different matrices must have the same name.
#'
#' @return A list with two elements:
#' \itemize{
#' \item WAP : the matrix resulting from WAP
#' \item alpha : a numeric vector with the weight coefficients of the networks
#' }
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Create random vector of accuracy metrics
#' m <- runif(3, min = 0, max = 1);
#'
#' # Integrate networks using Weighted Average Per-class (WAP) method
#' A_int <- WAP.int(m, align=TRUE, logint=FALSE, A1, A2, A3);
WAP.int <- function(m, align=FALSE, logint=FALSE, ...) {
    n <- length(list(...))
    if (length(m)!=n)
        stop("WAP.int: number of matrices and number of the values of the metric do not match")
    if (align)
        nets <- align.networks(fill=0, ...)
    if (logint) {
        m[m>0.99] <- 0.99
        m <- -log(1 - m)
    }
    norm.factor <- sum(m)
    alpha <- m/norm.factor
    res <- nets[[1]] * alpha[1]
    for (i in 2:n)
        res <- res + (nets[[i]] * alpha[i])
    return(list(WAP=res, alpha=alpha))
}



#' Weighted Average (WA) network integration
#'
#' @description It performs the WA integration between networks.
#'
#' Note that this function puts more weight (alpha parameter) for networks with
#' associated larger M. The alphas are computed by averaging across the alpha
#' of each class, and hence a unique integrated network is available for all the
#' considered classes.
#'
#' @param M a numeric matrix with the values of the metric used to compute the alpha
#'    coefficients for each class and for each network.
#'    Rows correspond to networks and columns to classes. Element (i,j) of the
#'    matrix corresponds to the value of the metric (e.g. AUC) for the ith
#'    network and the jth class.
#' @param logint logic. If TRUE the mean values m are log transformed: -log(1-m),
#'    otherwise a linear integration is performed (def: FALSE).
#' @param ... a list of numeric matrices. These must be named matrices representing
#'    adjacency matrices of the networks. Matrices may have different dimensions,
#'    but corresponding elements in different matrices must have the same name.
#'
#' @return A list with two elements:
#' \itemize{
#' \item WA : the matrix resulting from WA
#' \item alpha : a numeric vector with the weight coefficients of the networks
#' }
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Create random matrix of accuracy metrics (considering 3 classes)
#' M <- matrix(runif(9, min = 0, max = 1), ncol = 3);
#'
#' # Integrate networks using Weighted Average (WA) method
#' A_int <- WA.int(M, logint=TRUE, A1, A2, A3);
WA.int <- function(M, logint=FALSE, ...) {
    n <- length(list(...))
    if (nrow(M)!=n)
        stop("WA.int: number of matrices and number of rows of M do not match")
    n.class <- ncol(M)
    nets <- align.networks(fill=0, ...)

    alpha <- apply(M, 1, function(x) {return(mean(x));})

    if (logint) {
        alpha[alpha>0.99] <- 0.99
        alpha <- -log(1 - alpha)
    }
    alpha <- alpha/sum(alpha)

    res <- nets[[1]] * alpha[1]
    for (i in 2:n)
        res <- res + (nets[[i]] * alpha[i])
    return(list(WA=res, alpha=alpha))
}



###############################################################################
################## Utilities ##################################################
###############################################################################

#' Function to align the adjacency matrices of graphs/networks
#'
#' @description It accepts a list of adjacency matrices (that is an arbitrary number
#' of matrices separated by commas) and returns another list of adjacency matrices
#' having as elements the union of the elements of all the matrices.
#' Missed elements are replaced with null rows and columns. In this way the
#' resulting matrices have the same number of rows/columns in the same order.
#'
#' @param fill value used for the missing elements (def: 0).
#' @param ... a list of numeric matrices. These must be named matrices, and
#'     corresponding elements in different matrices must have the same name.
#'
#' @return  list of matrices : they correspond exactly and in the same order to
#' the input matrices, but they are filled with rows and columns when they have
#' missed values. The missing values are filled with fill.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Align networks
#' A_aligned <- align.networks(fill = 0, A1, A2, A3);
align.networks <- function(fill=0, ...) {
    networks <- list(...)

    n <- length(networks)
    if (n < 2)
        stop("align.networks: more than 1 network matrix must be provided")

    elem <- character()
    for (i in 1:n) {
        names.elem <- rownames(networks[[i]])
        if (is.null(names.elem))
            stop("align.networks: network matrices must be named")
        elem <- union(elem,names.elem)
    }

    n.elem <- length(elem)
    list.net <- list()
    for (i in 1:n) {
        m <- networks[[i]]
        mm <- matrix(rep(fill,n.elem*n.elem), nrow=n.elem)
        rownames(mm) <- colnames(mm) <- elem
        mm[rownames(m),colnames(m)] <- m
        list.net[[i]] <- mm
    }
    return(list.net)
}



#' Function to align a list of adjacency matrices of graphs/networks
#'
#' @description It accepts a list of adjacency matrices (that is an arbitrary number
#' of matrices separated by commas) and returns another list of adjacency matrices
#' having as elements the union of the elements of all the matrices.
#' Missed elements are replaced with null rows and columns. In this way the
#' resulting matrices have the same number of rows/columns in the same order.
#' NOTE: It is equal to align.networks with a list of matrices as argument instead of
#' the ... generic argument.
#'
#' @param fill value used for the missing elements (def: 0).
#' @param networks a list of numeric matrices. These must be named matrices, and
#'     corresponding elements in different matrices must have the same name.
#'
#' @return list of matrices : they correspond exactly and in the same order to the input
#' matrices, but they are filled with rows and columns when they have missed values.
#' The missing values are filled with fill.
#' @export
#'
#' @examples
#' # Create three example networks of different size
#' set.seed(123);
#' A1 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A1[lower.tri(A1)] = t(A1)[lower.tri(A1)];
#' diag(A1) <- 0;
#' rownames(A1) <- colnames(A1) <- sample(LETTERS, 10);
#'
#' A2 <- matrix(runif(49, min = 0, max = 1), nrow = 7);
#' A2[lower.tri(A2)] = t(A2)[lower.tri(A2)];
#' diag(A2) <- 0;
#' rownames(A2) <- colnames(A2) <- rownames(A1)[1:7];
#'
#' A3 <- matrix(runif(100, min = 0, max = 1), nrow = 10);
#' A3[lower.tri(A3)] = t(A3)[lower.tri(A3)];
#' diag(A3) <- 0;
#' rownames(A3) <- colnames(A3) <- c(rownames(A1)[1:5], c("A", "B", "Z", "K", "Q"));
#'
#' # Align networks
#' A_aligned <- lalign.networks(fill = 0, list(A1, A2, A3));
lalign.networks <- function(fill=0, networks) {

    n <- length(networks)
    if (n < 2 )
        stop("align.networks: more than 1 network matrix must be provided")

    elem <- character()
    for (i in 1:n) {
        names.elem <- rownames(networks[[i]])
        if (is.null(names.elem))
            stop("align.networks: network matrices must be named")
        elem <- union(elem,names.elem)
    }

    n.elem <- length(elem)
    list.net <- list()
    for (i in 1:n) {
        m <- networks[[i]]
        mm <- matrix(rep(fill,n.elem*n.elem), nrow=n.elem)
        rownames(mm) <- colnames(mm) <- elem
        mm[rownames(m),colnames(m)] <- m
        list.net[[i]] <- mm
    }
    return(list.net)
}


