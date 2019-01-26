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
X3 <- data.frame(x1 = rnorm(500, mean = -8, sd = 1),
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
    maxx1 <- c()
    minx1 <- c()
    maxx2 <- c()
    minx2 <- c()
    for (i in (1:n_cla)) {
        maxx1[i] <- mu[i, 1] + 4 * sqrt(cov[1, 1])
        minx1[i] <- mu[i, 1] - 4 * sqrt(cov[1, 1])
        maxx2[i] <- mu[i, 2] + 4 * sqrt(cov[2, 2])
        minx2[i] <- mu[i, 2] - 4 * sqrt(cov[2, 2])
    }
    maxx1 <- max(maxx1)
    minx1 <- min(minx1)
    maxx2 <- max(maxx2)
    minx2 <- min(minx2)
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
    for (i in (1:n_cla)) {
        points(mu[i,1], mu[i,2], pch = 19)
        for (k in (1:n_cla)) {
            points(0.5 * (mu[i,1] + mu[k,1]), 0.5 * (mu[i,2] + mu[k,2]), pch = 19)
        }
    }
    return(list(pi, cov, out))
}
lda2(c(-2, 5), X)
