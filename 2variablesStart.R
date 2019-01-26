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
X4 <- data.frame(
    x1 = rnorm(200, mean = -2, sd = 1),
    x2 = rnorm(200, mean = 0, sd = 1),
    class = "red"
)
X <- rbind(X1, X2, X3, X4)

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
    col <- heat.colors(n = n_cla)
    plot(
        X$x1[X$class == nam[1]],
        X$x2[X$class == nam[1]],
        xlim = c(min(minx[, 1]), max(maxx[, 1])),
        ylim = c(min(minx[, 2]), max(maxx[, 2])),
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
            orth[[i + k - 1]] <- c() # does this indexing work for >3 classes (empty elements?)
            orth[[i + k - 1]] <- inv %*% (mu[i,] - mu[k,])
        }
    }
    orth[[1]] <- NULL
    orth[[0.5 * n_cla * (n_cla - 1) + 1]] <- NULL # we have n(n-1)/2 segments between n points and n(n-1)/2+2 vectors in orth
    slope <- c()
    for (i in (1:(0.5 * n_cla * (n_cla - 1)))) {
        slope[i] <- -orth[[i]][1, ] / orth[[i]][2, ]
    }
    interc <- matrix(nrow = n_cla, ncol = n_cla)
    for (i in (1:n_cla)) {
        for (k in (i:n_cla)) {
            interc[i,k] <- -(slope[i] * 0.5 * (mu[i, 1] + mu[k, 1]) - 0.5 * (mu[i, 2] + mu[k, 2])) # intercept wil have to be double indexed (matrix, because of combinatorics)
        }
    }
    abline(b = slope[2],
           a = interc[3,3]) # [1,2], [2,3], [3,2] are the correct intercepts for three classes
    return(list(pi, cov, out, slope, interc))
}
lda2(c(-2, 5), X)
abline(a = 7.2, b = 0.835)
