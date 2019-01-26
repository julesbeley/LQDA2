X1 <-
    data.frame(
        x1 = rnorm(100, mean = 0, sd = 1),
        x2 = rnorm(100, mean = 3, sd = 1),
        class = "black"
    )
X2 <-
    data.frame(
        x1 = rnorm(200, mean = -4, sd = 1),
        x2 = rnorm(200, mean = 8, sd = 1),
        class = "white"
    )
X3 <- data.frame(
    x1 = rnorm(500, mean = -8, sd = 1),
    x2 = rnorm(500, mean = 0, sd = 1),
    class = "orange"
)
X <- rbind(X1, X2, X3)

lda2 <- function(value, data) {
    if (is.numeric(data[, 1] & is.numeric(data[, 2]))) {
        names(data) <- c("x1", "x2", "class")
    }
    tab <- table(data$class)
    n_cla <- dim(tab)
    nam <- names(tab)
    pi <- c()
    mu <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        pi[i] <- tab[i] / dim(data)[1]
        mu[i, 1] <- mean(data$x1[data$class == nam[i]]) # mean matrix: (class1:x1 ... classn:x1)
        mu[i, 2] <- mean(data$x2[data$class == nam[i]]) #              (class1:x2 ... classn:x2)
    }
    names(pi) <- nam
    val <- list()
    dev <- list()
    mul <- list()
    ccov <- list()
    for (i in (1:n_cla)) {
        val[[i]] <- matrix(nrow = tab[i], ncol = 2)
        val[[i]] <- cbind(data$x1[data$class == nam[i]], data$x2[data$class == nam[i]])
        dev[[i]] <- matrix(nrow = tab[i], ncol = 2)
        for (j in (1:tab[i])) {
            dev[[i]][j, ] <- val[[i]][j,] - mu[i,]
        }
        mul[[i]] <- matrix(nrow = 2, ncol = 2)
        ccov[[i]] <- matrix(nrow = 2, ncol = 2)
        mul[[i]] <- t(dev[[i]]) %*% dev[[i]]
        ccov[[i]] <- mul[[i]] / (dim(data)[1] - n_cla)
    }
    dis <- c()
    exp <- c()
    cov <- Reduce('+', ccov)
    rownames(cov) <- c("x1", "x2")
    colnames(cov) <- c("x1", "x2")
    inv <- solve(cov)
    for (i in (1:n_cla)) {
        dis[i] <- t(value) %*% inv %*% mu[i,] - 0.5 %*% t(mu[i,]) %*% inv %*% mu[i,] + log(pi[i])
    }
    prob <- c()
    prob <- exp(dis) / sum(exp(dis))
    names(dis) <- nam
    class <- names(dis)[dis == max(dis)]
    score <- max(dis)
    proba <- max(prob)
    out <- as.data.frame(cbind(class, score, proba))
    maxx <- matrix(nrow = n_cla, ncol = 2)
    minx <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        for (j in (1:2)) {
            maxx[i, j] <- mu[i, j] + 4 * sqrt(cov[j, j])
            minx[i, j] <- mu[i, j] - 4 * sqrt(cov[j, j])
        }
    }
    maxx1 <- max(maxx[, 1])
    minx1 <- min(minx[, 1])
    maxx2 <- max(maxx[, 2])
    minx2 <- min(minx[, 2])
    col <- heat.colors(n = n_cla)
    plot(
        X$x1[X$class == nam[1]],
        X$x2[X$class == nam[1]],
        xlim = c(minx1, maxx1),
        ylim = c(minx2, maxx2),
        col = col[1],
        xlab = "X1",
        ylab = "X2"
    )
    for (i in (2:n_cla)) {
        points(X$x1[X$class == nam[i]],
               X$x2[X$class == nam[i]],
               col = col[i])
    }
    orth <- list()
    for (i in (1:n_cla)) {
        points(mu[i, 1], mu[i, 2], pch = 19)
        for (k in (1:n_cla)) {
            points(0.5 * (mu[i, 1] + mu[k, 1]), 0.5 * (mu[i, 2] + mu[k, 2]), pch = 19)
            orth[[i + k - 1]] <- c()
            orth[[i + k - 1]] <- inv %*% (mu[i,] - mu[k,])
        }
    }
    orth[[1]] <- NULL
    orth[[0.5 * n_cla * (n_cla - 1) + 1]] <- NULL # we have n(n-1)/2 segments between n points and n(n-1)/2+2 vectors in orth
    slope <- c()
    for (i in (1:(0.5 * n_cla * (n_cla - 1)))) {
        
    }
    slope <- -orth[[1]][1, ] / orth[[1]][2, ]
    abline(b = slope,
           a = -(slope * 0.5 * (mu[1, 1] + mu[2, 1]) - 0.5 * (mu[2, 2] + mu[1, 2])))
    return(list(pi, cov, out, mu))
}
lda2(c(-2, 5), X)
abline(a = 7.2, b = 0.835)
