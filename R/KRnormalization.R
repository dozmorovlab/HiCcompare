#' Performs KR (Knight-Ruiz) normalization on a Hi-C matrix
#'
#' @export
#' @param A the matrix to be normalized - any columns/rows of
#'     0's will be removed before being normalized.
#' @details Performs KR normalization. The function is a translation
#'     of the `Matlab` code provided in the 2012 manuscript.
#'     Knight PA, Ruiz D. A fast algorithm for matrix balancing. IMA
#'     Journal of Numerical Analysis. Oxford University Press; 2012; drs019.
#' @return A KR normalized matrix
#' @examples
#' m <- matrix(rpois(100, 5), 10, 10)
#' KRnorm(m)


KRnorm = function(A) {
  # remove any cols/rows of 0s
  zeros = unique(which(colSums(A) == 0), which(rowSums(A) == 0))
  if (length(zeros) > 0) {
    A = A[-zeros, -zeros]
    message(paste0('Cols/Rows removed: '))
    message(paste(" ", zeros, sep = " "))
  }
  # initialize
  tol = 1e-6; delta = 0.1; Delta = 3; fl = 0;
  # change NAs in matrix to 0's
  NAlist = which(is.na(A), arr.ind = TRUE)
  A[is.na(A)] = 0
  n = nrow(A)
  e = matrix(1, nrow=n, ncol = 1)
  x0 = e
  res = matrix(nrow = n, ncol=1)
  # inner stopping criterior
  g=0.9; etamax = 0.1;
  eta = etamax; stop_tol = tol*.5;
  x = x0; rt = tol^2; v = x*(A %*% x); rk = 1 - v;
  rho_km1 = t(rk) %*% rk; rout = rho_km1; rold = rout;
  MVP = 0; # We'll count matrix vector products.
  i = 0; # Outer iteration count.
  while(rout > rt) { # Outer iteration
    i = i + 1; k = 0; y = e;
    innertol = max(c(eta^2 %*% rout, rt));
    while( rho_km1 > innertol ) { #Inner iteration by CG
      k = k +1
      if( k ==1) {
        z = rk / v; p=z; rho_km1 = t(rk)%*%z;
      }else {
        beta = rho_km1 %*% solve(rho_km2)
        p = z + beta*p
      }

      # update search direction efficiently
      w = x * (A%*%(x*p)) + v*p
      alpha = rho_km1 %*% solve(t(p) %*% w)
      ap = c(alpha) * p
      # test distance to boundary of cone
      ynew = y + ap;
      if(min(ynew) <= delta) {
        if(delta == 0) break()
        ind = which(ap < 0);
        gamma = min((delta - y[ind])/ap[ind]);
        y = y + gamma %*% ap;
        break()
      }
      if(max(ynew) >= Delta) {
        ind = which(ynew > Delta);
        gamma = min((Delta-y[ind])/ap[ind]);
        y = y + gamma %*% ap;
        break()
      }
      y = ynew;
      rk = rk - c(alpha) * w; rho_km2 = rho_km1;
      Z = rk/v; rho_km1 = t(rk) %*% z;
    }
    x = x*y; v = x*(A %*% x);
    rk = 1 - v; rho_km1 = t(rk) %*% rk; rout = rho_km1;
    MVP = MVP + k + 1;
    # Update inner iteration stopping criterion.
    rat = rout %*% solve(rold); rold = rout; res_norm = sqrt(rout);
    eta_o = eta; eta = g %*% rat;
    if(g %*% eta_o^2 > 0.1) {
      eta = max(c(eta, g %*% eta_o^2));
    }
    eta = max(c(min(c(eta, etamax)), stop_tol/res_norm));
  }
  result = diag(c(x)) %*% A %*% diag(c(x))
  # reintroduce NAs in final matrix
  if(nrow(NAlist) > 0) {
    # for(i in 1:nrow(NAlist)) {
    #   result[NAlist[i,1], NAlist[i,2]] = NA
    # }
    idx <- as.matrix(NAlist[, 1:2])
    result[idx] <- NA
  }
  return(result)
}
