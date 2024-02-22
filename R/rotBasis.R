#' @title Generate a complete orthonormal basis including vectors orthogonal to a given vector
#' @description With its first column being the given vector divided by its norm, 
#' generate an orthogonal matrix by rotating the identity matrix. 
#' @param x A vector of length \eqn{d} whose direction will be the first column 
#' of the resulting orthogonal matrix. 
#' @details Let \eqn{e_k} be the \eqn{k}th column of a \eqn{d}-dimensional identity matrix, 
#' for \eqn{k=1,\dots,d}. 
#' With \eqn{e_x = ( x - (x^\top e_1) e_1 ) / \| x - (x^\top e_1) e_1 \|} such that 
#' \eqn{x = (x^\top e_1) e_1 + (x^\top e_x) e_x}, 
#' let \eqn{Q} be the Q component of the QR decomposition of 
#' \eqn{( e_1\ e_x\ e_2\ \cdots\ e_{d-1})}. 
#' Then the output orthogonal matrix is given by \eqn{Q D Q^\top}. 
#' Here, \eqn{D} is a block diagonal matrix with two diagonal blocks.
#' The first diagonal block is a \eqn{2}-by-\eqn{2} matrix with entries given by 
#' \eqn{D_{11} = x^\top e_1}, \eqn{D_{12} = -x^\top e_x}, 
#' \eqn{D_{21} = x^\top e_x}, and \eqn{D_{22} = x^\top e_1}. 
#' The second diagonal block is a \eqn{(d-2)}-dimensional identity matrix. 
#' @return A \eqn{d}-by-\eqn{d} orthogonal matrix, 
#' of which the first column being the direction of \code{x}, 
#' a complete orthonormal basis consists of columns of the output matrix.
#' @examples
#' rotBasis(c(1,0,0,0))
#' rotBasis(c(2,0,0,0))
#' rotBasis(c(1,-1,0,0))
#' @export
#' 
rotBasis <- function(x) {
  xnorm <- sqrt( sum(x^2) )
  if ( !isTRUE( all.equal( xnorm, 1 ) ) ) {
    x <- x / xnorm
  }
  d <- length(x)
  basisIn <- diag(d)
  Qmat <- qr.Q( qr( cbind( basisIn[,1], x, basisIn[,-1] ) ) )
  if ( ! isTRUE( all.equal( Qmat[,1], basisIn[,1] ) ) ) {
    Qmat <- -Qmat
  }
  angle <- c( x %*% Qmat[,1:2] )
  Rmid <- diag( c( rep(angle[1],2), rep(1,d-2) ) )
  Rmid[1,2] <- -angle[2]
  Rmid[2,1] <- angle[2]
  Rmat <- Qmat %*% Rmid %*% t(Qmat)
  return ( Rmat %*% basisIn )
}
