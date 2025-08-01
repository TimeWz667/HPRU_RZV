
find_lor <- function(p0, p1) log(p1 / (1 - p1)) - log(p0 / (1 - p0))


apply_lor <- function(p0, lor) {
  lo0 <- log(p0 / (1 - p0))
  1 / (1 + exp(- lo0 - lor))
}


align_to <- function(ves0, tar) {
  tar <- 0.9
  ves0 <- runif(300, 0.5, 0.7)
  lor <- find_lor(mean(ves0), tar)
  lor
  
  ves1 <- apply_lor(ves0, lor)
  return(ves1)
}


get_or <- function(mlu, ref_mlu) {
  
  or0 <- (mlu[1] / (1 - mlu[1])) / (ref_mlu[1] / (1 - ref_mlu[1]))
  
  n_run <- 1e5
  while (T) {
    p0 <- runif(1e5, ref_mlu[2], ref_mlu[3])
    lors <- rnorm(1e5, 0, 10)
    
    p1 <- apply_lor(p0, lors)
    sel <- (p1 < mlu[3]) & (p1 > mlu[2])
    
    if (mean(sel) > 0.02 | n_run > 1e7) {
      break
    } else {
      n_run <- round(n_run * 10)  
    }
  }
  return(c(or0, quantile(exp(lors[sel]), c(0.025, 0.975))))
}
