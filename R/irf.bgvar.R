irf.bgvar2 <- function (x, n.ahead = 24, shockinfo = NULL, quantiles = NULL, 
          expert = NULL, verbose = TRUE) 
{
  start.irf <- Sys.time()
  ident <- attr(shockinfo, "ident")
  if (is.null(ident)) {
    ident <- "chol"
  }
  if (!ident %in% c("chol", "girf", "sign")) {
    stop("Please choose available identification scheme!")
  }
  if (is.null(shockinfo) && ident == "sign") {
    stop("Please provide 'shockinfo' argument.")
  }
  if (is.null(quantiles)) {
    quantiles <- c(0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95)
  }
  if (!is.numeric(quantiles)) {
    stop("Please provide 'quantiles' as numeric vector.")
  }
  if (!is.null(shockinfo)) {
    shockinfo <- shockinfo[!duplicated(shockinfo), ]
  }
  if (verbose) 
    cat("Start computing impulse response functions of Bayesian Global Vector Autoregression.\n\n")
  lags <- x$args$lags
  pmax <- max(lags)
  xglobal <- x$xglobal
  Traw <- nrow(xglobal)
  bigK <- ncol(xglobal)
  bigT <- Traw - pmax
  A_large <- x$stacked.results$A_large
  F_large <- x$stacked.results$F_large
  S_large <- x$stacked.results$S_large
  Ginv_large <- x$stacked.results$Ginv_large
  F.eigen <- x$stacked.results$F.eigen
  thindraws <- length(F.eigen)
  Global <- FALSE
  if (!is.null(shockinfo)) 
    Global <- ifelse(any(shockinfo$global), TRUE, FALSE)
  Rmed <- NULL
  rot.nr <- NULL
  xdat <- xglobal[(pmax + 1):Traw, , drop = FALSE]
  varNames <- colnames(xdat)
  cN <- unique(sapply(strsplit(varNames, ".", fixed = TRUE), 
                      function(x) x[1]))
  vars <- unique(sapply(strsplit(varNames, ".", fixed = TRUE), 
                        function(x) x[2]))
  N <- length(cN)
  Q <- length(quantiles)
  expert.list <- list(MaxTries = 100, save.store = FALSE, use_R = FALSE, 
                      applyfun = NULL, cores = NULL)
  if (!is.null(expert)) {
    if (!(is.null(expert$cores) || is.numeric(expert$cores) || 
          expert$cores %in% c("all", "half"))) {
      stop("Please provide the expert argument 'cores' in appropriate form. Please recheck.")
    }
    for (n in names(expert)) expert.list[[n]] <- expert[[n]]
  }
  MaxTries <- expert.list$MaxTries
  save.store <- expert.list$save.store
  use_R <- expert.list$use_R
  applyfun <- expert.list$applyfun
  cores <- expert.list$cores
  if (ident == "chol") {
    if (verbose) {
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
    }
    if (is.null(shockinfo)) {
      shockinfo <- get_shockinfo("chol", nr_rows = length(varNames))
      shockinfo$shock <- varNames
    }
    if (!all(c("shock", "scale") %in% colnames(shockinfo))) {
      stop("Please provide appropriate dataframe for argument 'shockinfo'. Respecify.")
    }
    if (!all(shockinfo$shock %in% varNames)) {
      stop("Please provide shock of 'shockinfo' only to variables available in the dataset used for estimation. Respecify.")
    }
    irf.fun <- .irf.chol
    shock.nr <- nrow(shockinfo)
    shocks <- shocknames <- unique(shockinfo$shock)
    select_shocks <- NULL
    for (ss in 1:shock.nr) select_shocks <- c(select_shocks, 
                                              which(shocks[ss] == varNames))
    scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    shock.cN <- sapply(strsplit(shockinfo$shock, ".", fixed = TRUE), 
                       function(x) x[1])
    shock.var <- sapply(strsplit(shockinfo$shock, ".", fixed = TRUE), 
                        function(x) x[2])
    shock.idx <- list()
    for (cc in 1:N) shock.idx[[cc]] <- grep(cN[cc], varNames)
    shock.cidx <- cN %in% shock.cN
    if (Global) {
      if (length(unique(shock.var[shockinfo$global])) != 
          1) {
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr - (sum(shockinfo$global) - 1)
      scale.new <- rep(1, shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.", unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for (ss in 1:shock.nr) {
        if (shockinfo$global[tt] == TRUE) {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shocks[shockinfo$global]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt <- max(which(shockinfo$global)) + 1
        }
        else {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shocks[tt]
          scale.new[ss] <- scale[tt]
          tt <- tt + 1
        }
      }
      scale <- scale.new
    }
    shocklist = list(shock.idx = shock.idx, shock.cidx = shock.cidx, 
                     plag = pmax, MaxTries = MaxTries)
  }
  else if (ident == "girf") {
    if (verbose) {
      cat("Identification scheme: Generalized impulse responses.\n")
    }
    if (is.null(shockinfo)) {
      shockinfo <- get_shockinfo("girf", nr_rows = length(varNames))
      shockinfo$shock <- varNames
    }
    if (!all(c("shock", "scale") %in% colnames(shockinfo))) {
      stop("Please provide appropriate dataframe for argument 'shockinfo'. Respecify.")
    }
    if (!all(shockinfo$shock %in% varNames)) {
      stop("Please provide shock of 'shockinfo' only to variables available in the dataset used for estimation. Respecify.")
    }
    if (!is.null(shockinfo)) {
      shocks <- shocknames <- unique(shockinfo$shock)
      scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    }
    else {
      shocks <- shocknames <- varNames
      scale <- rep(1, length(shocks))
    }
    irf.fun <- .irf.girf
    shock.nr <- length(shocks)
    select_shocks <- NULL
    for (ss in 1:shock.nr) select_shocks <- c(select_shocks, 
                                              which(shocks[ss] == varNames))
    shock.idx <- list()
    for (cc in 1:N) shock.idx[[cc]] <- grep(cN[cc], varNames)
    shock.cidx <- rep(FALSE, N)
    if (Global) {
      shock.var <- sapply(strsplit(shockinfo$shock, ".", 
                                   fixed = TRUE), function(x) x[2])
      if (length(unique(shock.var[shockinfo$global])) != 
          1) {
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr - (sum(shockinfo$global) - 1)
      scale.new <- rep(1, shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.", unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for (ss in 1:shock.nr) {
        if (shockinfo$global[tt] == TRUE) {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shocks[shockinfo$global]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt <- max(which(shockinfo$global)) + 1
        }
        else {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shocks[tt]
          scale.new[ss] <- scale[tt]
          tt <- tt + 1
        }
      }
      scale <- scale.new
    }
    shocklist = list(shock.idx = shock.idx, shock.cidx = shock.cidx, 
                     plag = pmax, MaxTries = MaxTries)
  }
  else if (ident == "sign") {
    if (!all(c("shock", "restriction", "sign", "horizon", 
               "scale", "prob") %in% colnames(shockinfo))) {
      stop("Please provide columns 'shock', 'restriction', 'sign', 'horizon' and 'scal' in dataframe 'shockinfo'.")
    }
    if (!(all(shockinfo$shock %in% varNames) && all(shockinfo$restriction %in% 
                                                    varNames))) {
      stop("Please provide in columns 'shock' and 'restriction' of 'shockinfo' only variable names available in the dataset used for estimation. Respecify.")
    }
    if (!any(shockinfo$sign %in% c(">", "<", "0", "ratio.H", 
                                   "ratio.avg"))) {
      stop("Misspecification in 'sign'. Only the following is allowed: <, >, 0, ratio.H, ratio.avg")
    }
    if (verbose) {
      cat("Identification scheme: identification via sign-restriction.\n")
    }
    irf.fun <- .irf.sign.zero
    shocks <- shocknames <- unique(shockinfo$shock)
    shock.nr <- length(shocks)
    select_shocks <- NULL
    for (ss in 1:shock.nr) select_shocks <- c(select_shocks, 
                                              which(shocks[ss] == varNames))
    shock.cN <- unique(sapply(strsplit(shockinfo$shock, ".", 
                                       fixed = TRUE), function(x) x[1]))
    shock.var <- sapply(strsplit(shockinfo$shock, ".", fixed = TRUE), 
                        function(x) x[2])
    shock.idx <- list()
    for (cc in 1:N) shock.idx[[cc]] <- grep(cN[cc], varNames)
    shock.cidx <- cN %in% shock.cN
    scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    if (Global) {
      if (length(unique(shock.var[shockinfo$global])) != 
          1) {
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr - (sum(shockinfo$global[!duplicated(shockinfo$shock)]) - 
                                1)
      scale.new <- rep(1, shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.", unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for (ss in 1:shock.nr) {
        if (shockinfo$global[tt] == TRUE) {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shocks[shockinfo$global[!duplicated(shockinfo$shock)]]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt <- max(which(shockinfo$global)) + 1
        }
        else {
          shock.global[[shocknames[ss]]] <- varNames %in% 
            shockinfo$shock[tt]
          scale.new[ss] <- scale[tt]
          tt <- tt + 1
        }
      }
      scale <- scale.new
    }
    if (any(shockinfo$sign %in% c("0", "ratio.H", "ratio.avg"))) {
      for (ss in 1:length(shocks)) {
        idx <- shockinfo$sign[grep(shocks[ss], shockinfo$shock)] %in% 
          c("0", "ratio.H", "ratio.avg")
        # if (!all(sapply(strsplit(shockinfo$restriction[grep(shocks[ss], 
        #                                                     shockinfo$shock)[idx]], ".", fixed = TRUE), 
        #                 function(x) x[1]) == shock.cN[ss])) {
        #   stop("Please provide zero and rationality conditions only in same country as the origin of the shock.")
        # }
      }
    }
    if (any(shockinfo$sign == "ratio.H")) {
      idx <- which(shockinfo$sign == "ratio.H")
      for (ii in idx) {
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock + 1):(Kshock + 2), ] <- NA
        shockinfo$shock[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$shock[ii], 
                                                             2)
        shockinfo$restriction[(Kshock + 1):nrow(shockinfo)] <- c(shockinfo$restriction[ii], 
                                                                 strsplit(shockinfo$restriction[ii], "_")[[1]][1])
        shockinfo$sign[(Kshock + 1):nrow(shockinfo)] <- c("0", 
                                                          "-1")
        shockinfo$horizon[(Kshock + 1):nrow(shockinfo)] <- c(1, 
                                                             Mshock)
        shockinfo$scale[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$scale[ii], 
                                                             2)
        shockinfo$prob[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$prob[ii], 
                                                            2)
      }
      shockinfo <- shockinfo[-idx, ]
      rownames(shockinfo) <- seq(1, nrow(shockinfo))
    }
    if (any(shockinfo$sign == "ratio.avg")) {
      idx <- which(shockinfo$sign == "ratio.avg")
      for (ii in idx) {
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock + 1):(Kshock + Mshock), ] <- NA
        shockinfo$shock[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$shock[ii], 
                                                             Mshock)
        shockinfo$restriction[(Kshock + 1):nrow(shockinfo)] <- c(shockinfo$restriction[ii], 
                                                                 rep(strsplit(shockinfo$restriction[ii], "_")[[1]][1], 
                                                                     Mshock - 1))
        shockinfo$sign[(Kshock + 1):nrow(shockinfo)] <- c("0", 
                                                          rep(-1/(Mshock - 1), Mshock - 1))
        shockinfo$horizon[(Kshock + 1):nrow(shockinfo)] <- seq(1, 
                                                               Mshock)
        shockinfo$scale[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$scale[ii], 
                                                             Mshock)
        shockinfo$prob[(Kshock + 1):nrow(shockinfo)] <- rep(shockinfo$prob[ii], 
                                                            Mshock)
      }
      shockinfo <- shockinfo[-idx, ]
      rownames(shockinfo) <- seq(1, nrow(shockinfo))
    }
    sign.horizon <- unique(shockinfo$horizon)
    sign.horizon <- sort(sign.horizon, decreasing = FALSE)
    sign.shockvars <- unique(shockinfo$shock)
    H.restr <- length(sign.horizon)
    N.restr <- bigK * H.restr
    S.cube <- matrix(0, N.restr, bigK)
    P.cube <- matrix(0, N.restr, bigK)
    Z.cube <- array(NA, c(N.restr, N.restr, bigK))
    dimnames(S.cube)[[1]] <- dimnames(Z.cube)[[1]] <- dimnames(Z.cube)[[2]] <- dimnames(P.cube)[[1]] <- paste(rep(varNames, 
                                                                                                                  H.restr), ".", rep(sign.horizon, each = bigK), sep = "")
    dimnames(S.cube)[[2]] <- dimnames(Z.cube)[[3]] <- dimnames(P.cube)[[2]] <- varNames
    for (vv in 1:length(varNames)) {
      Z.temp <- matrix(0, N.restr, N.restr)
      if (varNames[vv] %in% sign.shockvars) {
        idx <- which(shockinfo$shock == varNames[vv])
        sign.restr <- shockinfo$restriction[idx]
        sign.signs <- shockinfo$sign[idx]
        sign.horiz <- shockinfo$horizon[idx]
        sign.probs <- shockinfo$prob[idx]
        s.point <- which(sign.signs == "<" | sign.signs == 
                           ">")
        z.point <- seq(1, length(idx))[-s.point]
        S.cube[paste(varNames[vv], ".1", sep = ""), varNames[vv]] <- 1
        P.cube[paste(varNames[vv], ".1", sep = ""), varNames[vv]] <- 1
        if (length(s.point) > 0) {
          for (ss in 1:length(s.point)) {
            S.cube[paste(sign.restr[s.point[ss]], sign.horiz[s.point[ss]], 
                         sep = "."), varNames[vv]] <- ifelse(sign.signs[s.point[ss]] == 
                                                               "<", -1, 1)
            P.cube[paste(sign.restr[s.point[ss]], sign.horiz[s.point[ss]], 
                         sep = "."), varNames[vv]] <- sign.probs[s.point[ss]]
          }
        }
        if (length(z.point) > 0) {
          for (zz in 1:length(z.point)) {
            if (sign.signs[z.point[zz]] == "0") {
              grp <- which(sign.horiz[z.point[zz]] == 
                             sign.horizon)
              row <- (grp - 1) * bigK + which(sign.restr[z.point[zz]] == 
                                                varNames)
              Z.temp[row, row] <- 1
            }
            else {
              grp <- which(sign.horiz[z.point[zz]] == 
                             sign.horizon)
              col <- (grp - 1) * bigK + which(sign.restr[z.point[zz]] == 
                                                varNames)
              Z.temp[row, col] <- as.numeric(sign.signs[z.point[zz]])
            }
          }
        }
      }
      Z.cube[, , vv] <- Z.temp
    }
    no.zero.restr <- rep(TRUE, N)
    shock.order <- seq(bigK)
    for (cc in 1:N) {
      idx <- shock.idx[[cc]]
      no.zero.restr[cc] <- ifelse(base::sum(abs(Z.cube[, 
                                                       , idx])) > 0, FALSE, TRUE)
      shock.names <- names(sort(apply(Z.cube[, , idx], 
                                      3, function(x) base::sum(abs(x))), decreasing = TRUE))
      for (kk in 1:length(shock.names)) shock.order[idx[kk]] <- which(shock.names[kk] == 
                                                                        varNames)
    }
    shocklist <- list(shock.idx = shock.idx, shock.cidx = shock.cidx, 
                      MaxTries = MaxTries, S.cube = S.cube, Z.cube = Z.cube, 
                      P.cube = P.cube, shock.order = shock.order, shock.horz = sign.horizon, 
                      plag = pmax, no.zero.restr = no.zero.restr)
    rm(S.cube, Z.cube, P.cube)
  }
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    }
    else {
      if (.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, 
                                                  X, FUN, ...)
      }
      else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, 
                                                 ..., mc.cores = cores)
      }
    }
  }
  if (is.null(cores)) 
    cores <- 1
  if (ident == "sign") {
    R_store <- array(NA_real_, dim = c(bigK, bigK, thindraws), 
                     dimnames = list(colnames(xglobal), colnames(xglobal), 
                                     NULL))
  }
  else {
    R_store <- NULL
  }
  IRF_store <- array(NA_real_, dim = c(bigK, bigK, n.ahead + 
                                         1, thindraws), dimnames = list(colnames(xglobal), paste0("shock_", 
                                                                                                  colnames(xglobal)), seq(0, n.ahead), NULL))
  imp_posterior <- array(NA_real_, dim = c(bigK, n.ahead + 
                                             1, shock.nr, Q))
  dimnames(imp_posterior) <- list(colnames(xglobal), seq(0, 
                                                         n.ahead), shocknames, paste0("Q", quantiles * 100))
  start.comp <- Sys.time()
  if (verbose) 
    cat(paste("Start impulse response analysis on ", cores, 
              " core", ifelse(cores > 1, "s", ""), " (", thindraws, 
              " stable draws in total).", sep = ""), "\n")
  if (use_R) {
    counter <- numeric(length = thindraws)
    imp.obj <- applyfun(1:thindraws, function(irep) {
      Ginv <- Ginv_large[, , irep]
      Fmat <- adrop(F_large[, , , irep, drop = FALSE], 
                    drop = 4)
      Smat <- S_large[, , irep]
      imp.obj <- irf.fun(xdat = xdat, plag = pmax, n.ahead = n.ahead, 
                         Ginv = Ginv, Fmat = Fmat, Smat = Smat, shocklist = shocklist)
      if (verbose) {
        if (ident == "sign") {
          if (!any(is.null(imp.obj$rot))) {
            cat("\n", as.character(Sys.time()), "MCMC draw", 
                irep, ": rotation found after ", imp.obj$icounter, 
                " tries", "\n")
          }
          else {
            cat("\n", as.character(Sys.time()), "MCMC draw", 
                irep, ": no rotation found", "\n")
          }
        }
      }
      return(list(impl = imp.obj$impl, rot = imp.obj$rot, 
                  icounter = imp.obj$icounter))
    })
    for (irep in 1:thindraws) {
      counter[irep] <- imp.obj[[1]]$icounter
      if (imp.obj[[1]]$icounter == MaxTries) {
        imp.obj[[1]] <- NULL
      }
      else {
        IRF_store[, , , irep] <- imp.obj[[1]]$impl
        if (ident == "sign") 
          R_store[, , irep] <- imp.obj[[1]]$rot
        imp.obj[[1]] <- NULL
      }
      if (irep%%50 == 0) 
        gc()
    }
    rm(imp.obj)
  }
  else {
    shocklist$shock.idx <- lapply(shocklist$shock.idx, function(l) l - 
                                    1)
    shocklist$shock.horz <- shocklist$shock.horz - 1
    shocklist$shock.order <- shocklist$shock.order - 1
    type <- ifelse(ident == "chol", 1, ifelse(ident == "girf", 
                                              2, 3))
    counter <- numeric(length = thindraws)
    save_rot <- ifelse(ident == "sign", TRUE, FALSE)
    temp = compute_irf(A_large = A_large, S_large = S_large, 
                       Ginv_large = Ginv_large, type = type, nhor = n.ahead + 
                         1, thindraws = thindraws, shocklist_in = shocklist, 
                       save_rot = save_rot, verbose = verbose)
    for (irep in 1:thindraws) {
      counter[irep] <- temp$counter[irep, 1]
      if (temp$counter[irep, 1] == MaxTries) {
        temp$irf[[1]] <- NULL
        temp$rot[[1]] <- NULL
      }
      else {
        IRF_store[, , , irep] <- temp$irf[[1]]
        if (ident == "sign") 
          R_store[, , irep] <- temp$rot[[1]]
        temp$irf[[1]] <- NULL
        temp$rot[[1]] <- NULL
      }
      if (irep%%50 == 0) 
        gc()
    }
    rm(temp)
    shocklist$shock.idx = lapply(shocklist$shock.idx, function(l) l + 
                                   1)
    shocklist$shock.horz = shocklist$shock.horz + 1
    shocklist$shock.order = shocklist$shock.order + 1
  }
  end.comp <- Sys.time()
  diff.comp <- difftime(end.comp, start.comp, units = "mins")
  mins <- round(diff.comp, 0)
  secs <- round((diff.comp - floor(diff.comp)) * 60, 0)
  if (verbose) 
    cat(paste("\nImpulse response analysis took ", mins, 
              " ", ifelse(mins == 1, "min", "mins"), " ", secs, 
              " ", ifelse(secs == 1, "second.\n", "seconds.\n"), 
              sep = ""))
  if (ident == "sign") {
    idx <- which(counter != MaxTries)
    rot.nr <- paste("For ", length(idx), " draws out of ", 
                    thindraws, " draws, a rotation matrix has been found.")
    if (length(idx) == 0) {
      stop("No rotation matrix found with imposed sign restrictions. Please respecify.")
    }
    if (verbose) 
      cat(rot.nr)
    Ginv_large <- Ginv_large[, , idx, drop = FALSE]
    A_large <- A_large[, , idx, drop = FALSE]
    S_large <- S_large[, , idx, drop = FALSE]
    thindraws <- length(idx)
  }
  if (Global) {
    impulse <- NULL
    for (ss in 1:shock.nr) {
      temp <- apply(IRF_store[, shock.global[[ss]], , , 
                              drop = FALSE], c(1, 3, 4), sum)
      Mean <- temp[which(shock.global[[ss]])[1], 1, ]
      for (irep in 1:thindraws) {
        temp[, , irep] <- (temp[, , irep]/Mean[irep]) * 
          scale[ss]
      }
      impulse <- abind(impulse, temp, along = 4)
    }
    IRF_store <- aperm(impulse, c(1, 4, 2, 3))
    dimnames(IRF_store)[[2]] <- names(shock.global)
  }
  else {
    IRF_store <- IRF_store[, select_shocks, , , drop = FALSE]
    for (ss in 1:shock.nr) {
      Mean <- IRF_store[select_shocks[ss], ss, 1, ]
      for (irep in 1:thindraws) {
        IRF_store[, ss, , irep] <- (IRF_store[, ss, , 
                                              irep]/Mean[irep]) * scale[ss]
      }
    }
  }
  for (ss in 1:shock.nr) {
    for (qq in 1:Q) {
      imp_posterior[, , ss, qq] <- apply(IRF_store[, ss, 
                                                   , ], c(1, 2), quantile, quantiles[qq], na.rm = TRUE)
    }
  }
  A <- apply(A_large, c(1, 2), median)
  Fmat <- apply(F_large, c(1, 2, 3), median)
  Ginv <- apply(Ginv_large, c(1, 2), median)
  Smat <- apply(S_large, c(1, 2), median)
  Sigma_u <- Ginv %*% Smat %*% t(Ginv)
  if (ident == "sign") {
    imp.obj <- try(irf.fun(xdat = xdat, plag = pmax, n.ahead = n.ahead, 
                           Ginv = Ginv, Fmat = Fmat, Smat = Smat, shocklist = shocklist), 
                   silent = TRUE)
    if (!is(imp.obj, "try-error")) {
      Rmed <- imp.obj$rot
    }
    else {
      Rmed <- NULL
    }
  }
  struc.obj <- list(A = A, Fmat = Fmat, Ginv = Ginv, Smat = Smat, 
                    Rmed = Rmed)
  model.obj <- list(xglobal = xglobal, lags = lags)
  out <- structure(list(posterior = imp_posterior, ident = ident, 
                        shockinfo = shockinfo, rot.nr = rot.nr, struc.obj = struc.obj, 
                        model.obj = model.obj), class = "bgvar.irf")
  if (save.store) {
    out$IRF_store = IRF_store
    out$R_store = R_store
  }
  if (verbose) 
    cat(paste("\nSize of irf object: ", format(object.size(out), 
                                               unit = "MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf, start.irf, units = "mins")
  mins.irf <- round(diff.irf, 0)
  secs.irf <- round((diff.irf - floor(diff.irf)) * 60, 0)
  if (verbose) 
    cat(paste("\nNeeded time for impulse response analysis: ", 
              mins.irf, " ", ifelse(mins.irf == 1, "min", "mins"), 
              " ", secs.irf, " ", ifelse(secs.irf == 1, "second.", 
                                         "seconds.\n"), sep = ""))
  return(out)
}
