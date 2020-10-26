# The order lasso methods:
# for parameter tuning: plasso, slasso, alasso, talasso -> c.sel
# results: alasso, talasso, plasso, slasso



# -----------------------------------------------------
# Estimation
# -----------------------------------------------------
simulation.est = function(y, X, Lambda, sel) {
    # require("LasForecast")
    
    p = ncol(X)
    
    lambda.ada = Lambda$adalasso
    lambda.lasso = Lambda$lasso
    lambda.lasso.std = Lambda$lasso_std
    lambda.rep = Lambda$replasso
    
    # OLS
    coef.ols = lsfit(X, y, intercept = TRUE)$coef
    
    # OLS oracle
    coef.ora.temp = lsfit(X[, sel], y, intercept = TRUE)$coef
    coef.ora = rep(0, p + 1)
    coef.ora[c(TRUE, sel)] = coef.ora.temp
    
    # Lasso
    lasso = glmnet(
        X,
        y,
        lambda = lambda.lasso,
        intercept = TRUE,
        standardize = FALSE
    )
    a.lasso = as.numeric(lasso$a0)
    b.lasso = as.numeric(lasso$beta)
    coef.lasso = c(a.lasso, b.lasso)
    
    # Lasso_std
    lasso.std = glmnet(
        X,
        y,
        lambda = lambda.lasso.std,
        intercept = TRUE,
        standardize = TRUE
    )
    a.lasso.std = as.numeric(lasso.std$a0)
    b.lasso.std = as.numeric(lasso.std$beta)
    coef.lasso.std = c(a.lasso.std, b.lasso.std)
    
    # Adalasso
    alasso = adalasso(X, y, lambda = lambda.ada)
    a.ada = as.numeric(alasso$ahat)
    b.ada = as.numeric(alasso$bhat)
    coef.ada = c(a.ada, b.ada)
    
    # Replasso
    rlasso = replasso(X, y, b.ada, lambda.rep)
    a.rep = as.numeric(rlasso$ahat)
    b.rep = as.numeric(rlasso$bhat)
    coef.rep = c(a.rep, b.rep)
    
    coef.est = cbind(coef.ols,
                     coef.ora,
                     coef.lasso,
                     coef.lasso.std,
                     coef.ada,
                     coef.rep)
    result = split(coef.est, rep(1:ncol(coef.est), each = nrow(coef.est)))
    names(result) = c("ols",
                      "oracle",
                      "lasso",
                      "lasso_std",
                      "adalasso",
                      "replasso")
    
    return(result)
}

# -----------------------------------------------------
# Simulation pre-parameter-tuning
# -----------------------------------------------------
cross_validate <- function(X, y, c_seq, sel, k = 10) {
    # Cross-Validation function for simulation studies
    
    t = length(y)
    p = ncol(X)
    # number of methods used: lasso, lasso_std, adalasso, replasso
    m = 4
    
    seq.interval = split(1:t, ceiling(seq_along(1:t) / (t / k)))
    
    MSE = matrix(0, length(c_seq), m)
    
    for (i in 1:length(c_seq)) {
        c_i = c_seq[i]
        Lambda = as.list(c(rep(c_i / sqrt(t), 2), rep(c_i / (sqrt(t) * log(log(t))), 2)))
        names(Lambda) = c("lasso", "lasso_std", "adalasso", "replasso")
        
        for (j in 1:k) {
            y.est = y[-seq.interval[[j]]]
            X.est = X[-seq.interval[[j]],]
            
            yp = matrix(y[seq.interval[[j]]], length(seq.interval[[j]]), m)
            Xp = X[seq.interval[[j]],]
            
            Coef.temp = simulation.est(y.est, X.est, Lambda, sel)
            
            coef.lasso = Coef.temp$lasso
            coef.lasso.std = Coef.temp$lasso_std
            coef.ada = Coef.temp$adalasso
            coef.rep = Coef.temp$replasso
            
            Coef = cbind(coef.lasso, coef.lasso.std, coef.ada, coef.rep)
            
            mse.j = colMeans((yp - cbind(rep(
                1, length(seq.interval[[j]])
            ), Xp) %*% Coef) ^ 2)
            MSE[i,] = MSE[i,] + mse.j
            
        }
        
    }
    
    ind.sel = apply(MSE, 2, which.min)
    c.cv = c_seq[ind.sel]
    
    return(c.cv)
    
}

fix.lambda = function(t,
                      R,
                      a0,
                      b0,
                      dgp,
                      sigma_w_mat,
                      phi,
                      ur_index,
                      c_seq = seq(1e-4, 1e-2, length.out = 2000),
                      k = 5) {
    # Require k-fold CV function and dgp function
    # Require adalasso function and multi_assign function
    # Require SparseM, glmnet Rmosek, doParallel, foreach packages
    
    # Use cross-validation to fix the constant c_lambda
    # for a DGP and sample size t
    
    # Run R repliactions obatin a seq of c_lambda for each method
    # Return the median as the fixed c_lambda
    
    p <- length(b0)
    sel <- (b0 != 0)
    
    c_chosen = matrix(0, R, 4)
    colnames(c_chosen) = c("lasso", "lasso_std", "adalasso", "replasso")
    
    
    num_cores <- parallel::detectCores()
    cltr <- makeCluster(num_cores)
    registerDoParallel(cltr, cores = num_cores)
    
    c_chosen = foreach(
        r = 1:R,
        .combine = rbind,
        .packages = c("LasForecast", "SparseM"),
        .export = c("cross_validate",
                    "simulation.est",
                    "generate_innovation")
    ) %dorng% {
        # rewrite as a if-else structure
        
        D = dgp(t, c(a0, b0), sigma_w_mat, phi, ur_index = ur_index)
        y = D$y_est
        x = D$x_est
        
        c.cv = cross_validate(x, y, c_seq, sel, k = k)
        c.cv
    }
    stopCluster(cltr)
    
    c_seq = apply(c_chosen, 2, median)
    
    return(c_seq)
}

# -----------------------------------------------------
# Simulation main function
# -----------------------------------------------------
run_simulation <- function(dgp,
                           sigma_w_mat,
                           phi,
                           a0,
                           b0,
                           t_range,
                           R,
                           m,
                           ml,
                           lasso_names,
                           method_names,
                           c.sel,
                           ur_index) {
    p <- length(b0)
    sel <- (b0 != 0)
    
    num_cores <- parallel::detectCores()
    cltr <- makeCluster(num_cores)
    registerDoParallel(cltr, cores = num_cores)
    
    result <- foreach(
        n = t_range,
        .final = function(x) {
            setNames(x, names(t_range))
        },
        .packages = c("glmnet", "SparseM", "LasForecast"),
        .export = c("simulation.est", "generate_innovation"),
        .combine = list,
        .multicombine = TRUE
    ) %dorng% {
        # Result containers
        Y.True <- Y.Pred <- matrix(0, R, m + ml)
        Coef.est <- array(0, c(R, p, ml))
        colnames(Y.True) <- colnames(Y.Pred) <- method_names
        dimnames(Coef.est)[[3]] <- lasso_names
        
        Lambda <- as.list(c.sel / c(rep(sqrt(n), 2), rep(sqrt(n) * log(log(n)), 2)))
        names(Lambda) <- c("lasso", "lasso_std", "adalasso", "replasso")
        for (r in 1:R) {
            # Generate data
            
            D <- dgp(n, c(a0, b0), sigma_w_mat, phi, ur_index = ur_index)
            y <- D$y_est
            X <- D$x_est
            yp <- D$y_pred
            Xp <- D$x_pred
            
            Y.True[r,] <- yp
            
            Coef.temp <- simulation.est(y, X, Lambda, sel)
            
            coef.ols <- Coef.temp$ols
            coef.ora <- Coef.temp$oracle
            coef.lasso <- Coef.temp$lasso
            coef.lasso.std <- Coef.temp$lasso_std
            coef.ada <- Coef.temp$adalasso
            coef.rep <- Coef.temp$replasso
            
            coef.lasso.post <- post_lasso(X, y, coef.lasso)
            coef.lasso.std.post <- post_lasso(X, y, coef.lasso.std)
            coef.ada.post <- post_lasso(X, y, coef.ada)
            coef.rep.post <- post_lasso(X, y, coef.rep)
            
            # Generate Forecast
            Y.Pred[r,] <- c(
                crossprod(c(1, Xp), coef.ora),
                crossprod(c(1, Xp), coef.ols),
                crossprod(c(1, Xp), coef.ada),
                crossprod(c(1, Xp), coef.rep),
                crossprod(c(1, Xp), coef.lasso),
                crossprod(c(1, Xp), coef.lasso.std),
                crossprod(c(1, Xp), coef.ada.post),
                crossprod(c(1, Xp), coef.rep.post),
                crossprod(c(1, Xp), coef.lasso.post),
                crossprod(c(1, Xp), coef.lasso.std.post)
            )
            # Combine estimated coefficient
            Coef.est[r, ,] <- cbind(coef.ada[-1],
                                    coef.rep[-1],
                                    coef.lasso[-1],
                                    coef.lasso.std[-1])
        }
        
        # Compute FMSE
        FMSE <- apply((Y.True - Y.Pred) ^ 2, 2, mean)
        
        # Coefficient analysis
        Sel.est <- (Coef.est != 0)
        Sel.True <-
            array(matrix(sel, R, p, byrow = TRUE), c(R, p, ml))
        Correct.ratio <- apply(Sel.True == Sel.est, 3, mean)
        Select.ratio <-
            1 - colMeans(apply((Sel.True - Sel.est) == 1, c(1, 3), sum) / sum(sel))
        Screen.ratio <-
            1 - colMeans(apply((Sel.True - Sel.est) == -1, c(1, 3), sum) / (p - sum(sel)))
        
        coef.analysis <- list(correct = Correct.ratio,
                              select = Select.ratio,
                              screen = Screen.ratio)
        
        list(
            y0 = Y.True,
            yhat = Y.Pred,
            coef = Coef.est,
            FMSE = FMSE,
            coef.analysis = coef.analysis
        )
    }
    stopCluster(cltr)
    
    return(result)
    
}

sum_result <- function(Result){
    tt = length(Result)
    FMSE = NULL
    COEF = NULL
    for(t in 1:tt){
        
        FMSE = rbind(FMSE, Result[[t]]$FMSE)
        coef.t = c(
            Result[[t]]$coef.analysis$correct,
            Result[[t]]$coef.analysis$select,
            Result[[t]]$coef.analysis$screen
        )
        COEF = rbind(COEF, coef.t)
    }
    return(list(fmse = FMSE,coef = COEF))
}


sum_coint_coef <- function(result, inactive_index) {
    m <- length(result)
    both_sel <- NULL
    one_sel <- NULL
    no_sel <- NULL
    
    for (n in 1:m) {
        result_n <- result[[n]]
        
        sel_n <- (result_n$coef[, inactive_index,] != 0)
        both_sel <-
            rbind(both_sel, colMeans((sel_n[, 1, ] + sel_n[, 2, ]) == 2))
        one_sel <-
            rbind(one_sel, colMeans((sel_n[, 1, ] + sel_n[, 2, ]) == 1))
        no_sel <-
            rbind(no_sel, colMeans((sel_n[, 1, ] + sel_n[, 2, ]) == 0))
    }
    return(cbind(no_sel, one_sel, both_sel))
}