XY_linear_regression <- function(x, y, e, predict.X=NULL) {
    mx <- mean(x)
    my <- mean(y)
    n <- length(x)
    sx <- sd(x)*sqrt((n-1)/n)
    sy <- sd(y)*sqrt((n-1)/n)
    
    cov <- sum(x*y)/n - mx*my
    corr <- cov/(sx*sy)
    print(paste("Correlation", corr))
    
    T <- corr*sqrt((n-1)/(1-corr^2))
    T.crit <- qt(1-e/2, n-2)
    print(paste("T:", T, " T.crit:", T.crit))
    if (abs(T) < T.crit) {
        print("The correlation is considered 0!")
    }
    
    a <- corr*(sy/sx)
    b <- my - a*mx
    print(paste("a:", a, " b:", b))
    
    yp <- a*x+b
    SST <- sum((y-my)^2)
    SSE <- sum((y-yp)^2)
    
    R.2 <- 1 - SSE/SST
    print(paste("Coefficient of determination", R.2))
    
    SXX <- sum((x-mx)^2)
    SYY <- sum((y-my)^2)
    SXY <- sum((x-mx)*(y-my))
    
    verr <- SSE/(n-2)
    print(paste("Error's variance: ", verr))
    va <- verr/SXX
    vb <- verr*sum(x^2)/(n*SXX)
    sa <- sqrt(va)
    sb <- sqrt(vb)
    
    a.low <- a - sa*T.crit
    a.high <- a + sa*T.crit
    b.low <- b - sb*T.crit
    b.high <- b + sb*T.crit
    print(paste("a confidence interval: ", a.low, " - ", a.high))
    print(paste("b confidence interval: ", b.low, " - ", b.high))
    
    # assure that a cannot be zero
    if (a.low < 0 & 0 < a.high) {
        F <- R.2*(n-2)/(1-R.2)
        F.crit <- qf(1-e/2, 1, n-2)
        print(paste("F:", F, " F.crit:", F.crit))
        if (F < F.crit) {
            print("a is considered 0!")
        }
    }
    
    if (!is.null(predict.X)) {
        sy <- sqrt(verr*(1+1/n+(predict.X-mx)^2/SXX))
        Y.center <- a * predict.X + b
        
        Y.low <- Y.center - sy * T.crit
        Y.high <- Y.center + sy * T.crit
        print(paste("Y confidence interval :", Y.low, "-", Y.center, "-", Y.high))
    }
}











