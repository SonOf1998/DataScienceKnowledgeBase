one_sample_z_test <- function(x, e, sd, m0, two.tail=TRUE) {
    mx <- mean(x)
    n <- length(x)
    
    T <- (mx-m0)*sqrt(n)/sd
    T.crit <- if (two.tail) {
        qnorm(1-e/2)
    } else {
        qnorm(1-e)
    }
    
    print(paste("T:", T, "T.crit:", T.crit))
    if (abs(T) < T.crit) {
        print("The proposed mean can be accepted")
    } else {
        print("The proposed mean has to be rejected")
    }
}

two_sample_z_test <- function(x, y, e, sd, two.tail=TRUE) {
    mx <- mean(x)
    my <- mean(y)
    n <- length(x)
    m <- length(y)
    
    T <- (mx-my)/sqrt(sd^2/n + sd^2/m)
    T.crit <- if (two.tail) {
        qnorm(1-e/2)
    } else {
        qnorm(1-e)
    }
    
    print(paste("T:", T, "T.crit:", T.crit))
    if (abs(T) < T.crit) {
        print("The means can be considered equal")
    } else {
        print("The means are not equal")
    }
}

one_sample_t_test <- function(x, e, m0, two.tail=TRUE) {
    mx <- mean(x)
    n <- length(x)
    s <- sd(x)
    print(s)
    T <- (mx-m0)*sqrt(n)/s
    T.crit <- if (two.tail) {
        qt(1-e/2, n-1)
    } else {
        qt(1-e, n-1)
    }
    
    print(paste("T:", T, "T.crit:", T.crit))
    if (abs(T) < T.crit) {
        print("The proposed mean can be accepted")
    } else {
        print("The proposed mean has to be rejected")
    }
}


two_sample_t_test <- function(x, y, e, two.tail=TRUE) {
    mx <- mean(x)
    my <- mean(y)
    vx <- var(x)
    vy <- var(y)
    n <- length(x)
    m <- length(y)
    print(paste("VX: ",vx, "VY: ",vy, " These should be reviewed by an F test!"))
    
    T <- (mx - my)*sqrt(n*m*(n+m-2)/(n+m))/sqrt((n-1)*vx + (m-1)*vy)
    T.crit <- if (two.tail) {
        qt(1-e/2, n+m-2)
    } else {
        qt(1-e, n+m-2)
    }
    
    print(paste("T:", T, "T.crit:", T.crit))
    if (abs(T) < T.crit) {
        print("The means can be considered equal")
    } else {
        print("The means are not equal")
    }
}


f_test <- function(x, y, e) {
    n <- length(x)
    m <- length(y)
    vx <- var(x)
    vy <- var(y)
    
    f <- NULL
    f.crit <- NULL
    
    if (vx > vy) {
        f <- vx/vy
        f.crit <- qf(1-e/2, n-1, m-1)
    } else {
        f <- vy/vx
        f.crit <- qf(1-e/2, m-1, n-1)
    }
    
    print(paste("f:", f, "f.crit:", f.crit))
    if (f < f.crit) {
        print("Variances can be considered equal")
    } else {
        print("Variances are definitely different")
    }
}


welch_test <- function(x, y, e, two.tail=TRUE) {
    mx <- mean(x)
    my <- mean(y)
    n <- length(x)
    m <- length(y)
    vx <- var(x)*(n-1)/n
    vy <- var(y)*(n-1)/n
    
    W <- (mx - my)/sqrt(vx/n + vy/m)
    c <- (vx/m) / (vx/n + vy/m)
    f.inv <- c^2/(m-1)+(1-c)^2/(n-1)
    f <- 1/f.inv
    W.crit <- if (two.tail) {
        qt(1-e/2, f)
    } else {
        qt(1-e, f)
    }
    
    print(paste("W:", W, "W.crit:", W.crit))
    if (abs(W) < W.crit) {
        print("The means can be considered equal")
    } else {
        print("The means are not equal")
    }
}

library("BoutrosLab.plotting.general")
one_sample_kolmogorov_smirnov_test <- function(x, e, mu, sigma) {
    x <- sort(x)
    n <- length(x)
    empirical_cdf <- function(v) {
        if (v <= x[1]) {
            return(0)
        }
        if (v > x[n]) {
            return(1)
        }
        
        value <- 0
        for (i in 1:(n-1)){
            if (x[i] < v) {
                value <- value + 1/n
            }
            else {
                break
            }
        }
        
        value
    }
    
    diff.max <- 0
    for (elem in x) {
        # it is recommended to look for the maximum difference at the
        # discontinous points of the empirical cdf
        value <- elem + 0.00001
        diff <- abs(empirical_cdf(value) - pnorm(value, mean = mu, sd = sigma))
        diff.max <- max(diff.max, diff)
    }
    
    crit <- ks.test.critical.value(n, 1-e) 
    print(paste("Max difference:", diff.max, "critical value:", crit))
    if (diff.max < crit) {
        print(paste("The given sample follows the normal distribution of N(",mu,", ",sigma,").",sep=""))   
    } else {
        print(paste("The given sample does not follow the normal distribtuion of N(",mu,", ",sigma,").",sep=""))
    }
}

# can be used for paired sample test (med0 = 0) (x = X - Y)
one_sample_wilcoxon_test <- function(x, med0, e, two.tail=TRUE) {
    d <- x - med0
    d <- d[d!=0]
    s <- sign(d)
    a <- abs(d)
    r <- rank(a, ties.method="average")
    sr <- s*r
    n <- length(d)
    ranksum <- n*(n+1)/2
    
    R.plus <- sum(sr[sr>0])
    R.crit <- if (two.tail) {
        qsignrank(e/2, n) - 1
    } else {
        qsignrank(e, n) - 1
    }
    
    print(paste("R.plus:", R.plus, ", R.crit:", R.crit, ", ranksum:", ranksum))
    if (R.crit < R.plus  & R.plus < ranksum - R.crit) {
        print("The proposed median can be accepted")
    } else {
        print("The proposed median has to be rejected")
    }
}



#Also known as: two_sample_wilcoxon test
mann_whitney_u_test <- function(x, y, e, two.tail=TRUE, force.normal=FALSE) {
    n <- length(x)
    m <- length(y)
    z <- c(x, y)
    z.ordered <- sort(z)
    R.x <- sum(match(x, z.ordered))
    R.y <- sum(match(y, z.ordered))
    U.x <- R.x - n*(n+1)/2
    U.y <- R.y - m*(m+1)/2
    
    if (n + m < 20 & !force.normal)
    {
        #critical value
        c <- if (two.tail) {
            qwilcox(e/2, n, m) - 1
        } else {
            qwilcox(e, n, m) - 1
        }
        
        print(paste("U.x:", U.x, ", U.y:", U.y, ", critical value:", c))
        U <- min(U.x, U.y)
        if (c < U & U < U.x + U.y - c) {
            print("Equal median and identical distribution can be assumed")
        } else {
            print("Equal median and identical distribution cannot be assumed")
        }
    }
    else {
        mx <- n*(n+m+1)/2
        sx <- sqrt(n*m*(n+m+1)/12)
        u <- (R.x - mx)/sx
        
        c <- if (two.tail) {
            qnorm(1-e/2)
        } else {
            qnorm(1-e)
        }
        
        print(paste("R.x:", R.x, "mx:", mx, "sx:", sx, "u:", u, "u.crit", c))
        if (abs(u) < c) {
            print("Equal median and identical distribution can be assumed")
        }
        else {
            print("Equal median and identical distribution cannot be assumed")
        }
    }
}


chi_square_goodness_of_fit_test <- function(x, intervals, e, mu, sigma) {
    r <- length(x)
    n <- sum(x)
    
    p <- pnorm(intervals[2, ], mean=mu, sd=sigma) - pnorm(intervals[1, ], mean=mu, sd=sigma)
    
    T <- 0
    for (i in 1:r) {
        T <- T + (x[i] - n * p[i])^2/(n * p[i])
    }
    
    c.crit <- qchisq(1-e, r-1)
    print(paste("T:", T, "c.crit:", c.crit))
    if (T < c.crit) {
        print(paste("The given sample follows the normal distribution of N(",mu,", ",sigma,").",sep=""))   
    } else {
        print(paste("The given sample does not follow the normal distribtuion of N(",mu,", ",sigma,").",sep=""))
    }
}


chi_square_test_of_homogenity <- function(x, y, e) {
    r <- length(x)
    n <- sum(x)
    m <- sum(y)
    
    T <- 0
    for (i in 1:r) {
        T <- T + n*m*(x[i]/n - y[i]/m)^2/(x[i] + y[i])
    }
    
    c.crit <- qchisq(1-e, r-1)
    print(paste("T:", T, "c.crit:", c.crit))
    if (T < c.crit) {
        print("The samples have the same distribution")   
    } else {
        print("The samples follow a different distribution")
    }
}


chi_square_test_of_independencce <- function(m, e) {
    k <- dim(m)[1]
    l <- dim(m)[2]
    
    n <- sum(m)
    p.xy <- m/n
    x <- rowSums(m)
    y <- colSums(m)
    p.x <- x/n
    p.y <- y/n
    
    T <- 0
    for (i in 1:k) {
        for (j in 1:l) {
            T <- T + n * (p.xy[i, j]-p.x[i]*p.y[j])^2/(p.x[i]*p.y[j])
        }
    }
    
    c.crit <- qchisq(1-e, (k-1)*(l-1))
    print(paste("T:", T, "c.crit:", c.crit))
    if (T < c.crit) {
        print("The attributes can be considered independent")   
    } else {
        print("The attributes are dependent")
    }
}

ANOVA <- function(e, ...) {
    groups <- list(...)
    p <- length(groups)
    
    group.means <- sapply(groups, mean)
    ns <- sapply(groups, length)
    n <- length(unlist(groups))
    global.mean <- mean(unlist(groups))
    
    Q.k <- 0
    Q.b <- 0
    for (i in 1:p) {
        Q.k <- Q.k + ns[i]*(group.means[i]-global.mean)^2
        for (j in 1:ns[i]) {
            Q.b <- Q.b + (groups[[i]][j]-group.means[i])^2
        }
    }
    
    f <- (n-p)*Q.k/((p-1)*Q.b)   
    f.crit <- qf(1-e/2, p-1, n-p)
    print(paste("f:", f, "f.crit:", f.crit))
    if (f < f.crit) {
        print("The means can be considered equal")   
    } else {
        print("The means are not equal")
    }
}


library('SuppDists')
friedman_test <- function(m, e, force.chisq=FALSE) {
    n <- dim(m)[1]
    p <- dim(m)[2]
    
    m.ranked <- t(apply(m, 1, rank))
    Rs <- colSums(m.ranked)
    
    
    T <- 0
    for (i in 1:p) {
        T <- T + 12/(n*p*(p+1))*(Rs[i]^2)
    }
    T <- T - 3*n*(p+1)
    
    if (n < 20 & !force.chisq) {
        T.crit <- qFriedman(1-e, p, n)
        print(paste("T:", T, " T.crit", T.crit))
        if (T < T.crit) {
            print("The variables has the same median")
        }
        else {
            print("The medians of the samples are not indentical")
        }
    } else {
        T.crit <- qchisq(1-e, p-1)
        print(paste("T:", T, " T.crit", T.crit))
        if (T < T.crit) {
            print("The variables has the same median")
        }
        else {
            print("The medians of the samples are not indentical")
        }
    }
}

kruskal_wallis_test <- function(eps, m){
    n <- dim(m)[1]
    p <- dim(m)[2]
    N <- n*p
    
    m.ranked <- rank(m)
    dim(m.ranked) = c(n,p)
    Rs <- colSums(m.ranked)
    
    
    T <- 0
    for (i in 1:p) {
        T <- T + 12/(N*(N+1))*(Rs[i]^2/n)
    }
    T <- T - 3*(N+1)
    
    T.crit <- qchisq(1-e, p-1)
    print(paste("T:", T, " T.crit", T.crit))
    if (T < T.crit) {
        print("The variables has the same median")
    }
    else {
        print("The medians of the samples are not indentical")
    }
}


