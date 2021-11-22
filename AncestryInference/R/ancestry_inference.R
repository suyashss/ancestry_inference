#' Infer ancestry
#'
#' Infer ancestry from VCF file
#' 
#' @param vcffile A \code{character}. Path to an uncompressed vcf file.
#' @param frqfile A \code{character}. Path to a frequency file.
#' @param eps A \code{numeric}. Bounds minimum and maximum ancestries.
#' @param verbose A \code{logical}. If set to TRUE, display additional information.
#' 
#' @importFrom plyr join
#' @export
infer_ancestry <- function(vcffile, frqfile = NULL, eps = 1e-04, verbose = F) {
    frqtable <- read.table(frqfile, header = T)
    vcffull <- read.table(vcffile, header = F, stringsAsFactors = F)
    colnames(vcffull)[1:9] <- c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    vcf <- vcffull[, -1 * c(3, 6, 7, 8, 9)]

    ob <- join(frqtable, vcf, type = "left")
    poplist <- c("AFR", "ASN", "AMR", "SAN", "EUR")
    freq_subset <- 1 - as.matrix(ob[, poplist])
    nextracols <- length(poplist) + 4  # CHR,POS,REF,ALT
    cols <- (nextracols + 1):dim(ob)[2]
    ans <- matrix(numeric(), nrow = length(cols), ncol = length(poplist))
    ans2 <- matrix(numeric(), nrow = length(poplist), ncol = length(cols))
    colnames(ans) <- poplist
    for (ii in seq_along(cols)) {
        if (verbose) {
            print(paste("Sample", ii))
        }
        gtvec_text <- ob[, cols[ii]]
        gtvec_text[is.na(gtvec_text)] <- "0|0"

        gtvec <- sapply(gtvec_text, .gt_texttonum)
        qvec <- .findq(gtvec, freq_subset, eps = eps)
        ans[ii, ] <- qvec
        ans2[, ii] <- qvec
    }
    return(ans)
}

.vcf2gt <- function(x) {
    strsplit(x, split = ":")[[1]][1]
}

#' @importFrom data.table fread
infer_ancestry_dt <- function(vcffile, frqfile = NULL, eps = 1e-04, verbose = F) {
    frqtable <- fread(frqfile)
    vcffull <- data.table(read.table(vcffile, header = F, stringsAsFactors = F))
    colnames(vcffull)[1:9] <- c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    vcf <- vcffull[, -1 * c(3, 6, 7, 8, 9), with = F]
    vcffrq <- merge(frqtable, vcf, all.x = T, by = c("CHR", "POS", "REF", "ALT"))
    poplist <- c("AFR", "ASN", "AMR", "SAN", "EUR")
    freq_subset <- 1 - as.matrix(vcffrq[, poplist, with = F])
    nextracols <- length(poplist) + 4  # CHR,POS,REF,ALT
    gtcols <- names(vcffrq)[(nextracols+1):ncol(vcffrq)]
    vcffrq[, `:=`((gtcols), lapply(.SD, .gt_texttonum_vec)), .SDcols = gtcols]
    # Inference
    inf_anc <- vcffrq[, lapply(.SD, .findq, freq_subset, eps), .SDcols = gtcols]
    inf_anc[, population := poplist]
    manc <- melt(inf_anc, id.vars = "population", variable.name = "sample")
    danc <- dcast(manc, sample~population)
    danc[, InferredAncestry := names(.SD)[max.col(.SD)], .SDcols = poplist]
    return(ret)
}

.gt_texttonum <- function(gtstr) {
    x <- strsplit(gtstr, ":")[[1]][1]
    ans <- 0
    if(x=="0/1" || x=="0|1" || x=="1/0" || x=="1|0"){ # Slightly faster tham %in%
        ans <- 1
    } else if(x=="1|1" || x=="1/1"){
        ans <- 2
    } else {
        ans <- 0
    }
    return(ans)
}

#' @importFrom stringr str_count
.gt_texttonum_vec <- function(gtstr){
    gts <- sapply(strsplit(gtstr, split = ":"), "[[", 1)
    alts <- str_count(gts, "1")
    alts[is.na(alts)] <- 0
    return(alts)
}


.boundf <- function(f) {
    f1 <- f
    f1[f1 < eps] <- eps
    f1[f1 > 1 - eps] <- 1 - eps
    return(f1)
}

.ll <- function(x, g, fmat) {
    # fmat=JxK
    qvec <- x  #c(x,1-sum(x))
    qmat <- matrix(qvec, nrow = length(qvec), ncol = 1, byrow = F)  #Kx1
    fbar <- .boundf(c(fmat %*% qmat))
    -1 * sum(g * log(fbar) + (2 - g) * log(1 - fbar))
}

.gradll <- function(x, g, fmat) {
    # fmat=JxK
    qvec <- x  #c(x,1-sum(x))
    qmat <- matrix(qvec, nrow = length(qvec), ncol = 1, byrow = F)  # Kx1
    fbar <- .boundf(c(fmat %*% qmat))  #Jx1
    gmat <- matrix(g, ncol = length(qvec), nrow = length(g), byrow = F)  # JxK
    fbarmat <- matrix(fbar, nrow = length(fbar), ncol = length(qvec), byrow = F)
    gradvec <- c(colSums(fmat * gmat/fbarmat + (2 - gmat) * (1 - fmat)/(1 - fbarmat)))
    -1 * gradvec
}

.issumq1 <- function(x, g, fmat) {
    sum(x) - 1
}

.issumq1_g <- function(x, g, fmat) {
    rep(1, length(x))
}

.init_em <- function(g, fmat) {
    K = dim(fmat)[2]
    xold = rep(1/K, K)
    xnew = xold
    gmat = matrix(g, ncol = K, nrow = length(g), byrow = F)  # JxK
    for (i in 1:200) {
        qmat = matrix(xold, nrow = K, ncol = 1, byrow = F)  # Kx1
        fbar = .boundf(c(fmat %*% qmat))  #Jx1
        fbarmat = matrix(fbar, nrow = length(fbar), ncol = K, byrow = F)
        qmat2 = matrix(xold, nrow = dim(fmat)[1], ncol = K, byrow = T)
        amat = qmat2 * fmat/fbarmat
        bmat = qmat2 * (1 - fmat)/(1 - fbarmat)
        xnew = colSums(gmat * amat + (2 - gmat) * bmat)/(2 * length(fbar))
        xnew = xnew/sum(xnew)
        # cat('EM:',i,xold,'--->',xnew,'=',sum(xnew),'\n')
        xold = .boundf(xnew)
    }
    # cat('Initial guess is',xnew,'\n')
    .boundf(xnew)
}

#' @importFrom nloptr nloptr
.findq <- function(g, fmat, eps = NULL) {
    K <- dim(fmat)[2]
    local_opts <- list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-07)
    opts <- list(algorithm = "NLOPT_LD_SLSQP", print_level = 0, maxeval = 10000, xtol_rel = 1e-07, ftol_abs = 1e-10)
    x0 <- .init_em(g, fmat)
    solobj = nloptr(x0 = x0, eval_f = .ll, eval_grad_f = .gradll, eval_g_eq = .issumq1, eval_jac_g_eq = .issumq1_g,
        lb = rep(eps, K), ub = rep(1 - eps, K), g = g, fmat = fmat, opts = opts)
    bestq <- solobj$solution
    return(bestq)
}
