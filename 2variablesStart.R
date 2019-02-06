X1 <-
    data.frame(
        x1 = rnorm(100, mean = runif(1,-10,10), sd = 1),
        x2 = rnorm(100, mean = runif(1,-10,10), sd = 1),
        class = "black"
    )
X2 <-
    data.frame(
        x1 = rnorm(200, mean = runif(1,-10,10), sd = 1),
        x2 = rnorm(200, mean = runif(1,-10,10), sd = 1),
        class = "blue"
    )
X3 <- data.frame(
    x1 = rnorm(200, mean = runif(1,-10,10), sd = 1),
    x2 = rnorm(200, mean = runif(1,-10,10), sd = 1),
    class = "orange"
)
X4 <- data.frame(
    x1 = rnorm(500, mean = runif(1,-10,10), sd = 1),
    x2 = rnorm(500, mean = runif(1,-10,10), sd = 1),
    class = "red"
)
X5 <- data.frame(
    x1 = rnorm(600, mean = runif(1,-10,10), sd = 1),
    x2 = rnorm(600, mean = runif(1,-10,10), sd = 1),
    class = "white"
)
X6 <- data.frame(
    x1 = rnorm(800, mean = runif(1,-10,10), sd = 1),
    x2 = rnorm(800, mean = runif(1,-10,10), sd = 1),
    class = "green"
)
X <- rbind(X1, X2, X3, X4, X5, X6)

lda2 <- function(data) {
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
        mu[i, 1] <- mean(data$x1[data$class == nam[i]]) 
        mu[i, 2] <- mean(data$x2[data$class == nam[i]])
    }
    names(pi) <- nam
    val <- list()
    dev <- list()
    mul <- list()
    ccov <- list()
    for (i in (1:n_cla)) {
        val[[i]] <- matrix(nrow = tab[i], ncol = 2)
        val[[i]] <-
            cbind(data$x1[data$class == nam[i]], data$x2[data$class == nam[i]])
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
    maxx <- matrix(nrow = n_cla, ncol = 2)
    minx <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        for (j in (1:2)) {
            maxx[i, j] <- mu[i, j] + 4 * sqrt(cov[j, j])
            minx[i, j] <- mu[i, j] - 4 * sqrt(cov[j, j])
        }
    }
    col <- heat.colors(n = n_cla, alpha = 0.7)
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
    runifx1 <-
        runif(200000,
              min = 1.2 * min(minx[, 1]),
              max = 1.2 * max(maxx[, 1]))
    runifx2 <-
        runif(200000,
              min = 1.2 * min(minx[, 2]),
              max = 1.2 * max(maxx[, 2]))
    runif <- cbind(runifx1, runifx2)
    dismc <- matrix(nrow = 200000, ncol = n_cla)
    for (h in (1:200000)) {
        for (i in (1:n_cla)) {
            dismc[h, i] <-
                t(runif[h,]) %*% inv %*% mu[i,] - 0.5 %*% t(mu[i,]) %*% inv %*% mu[i,] + log(pi[i])
        }
    }
    colnames(dismc) <- nam
    rownames(mu) <- nam
    classmc <- c()
    for (h in (1:200000)) {
        classmc[h] <- names(dismc[h,])[dismc[h,] == max(dismc[h,])]
    }
    dismc <-
        data.frame(
            x1 = runifx1,
            x2 = runifx2,
            class = classmc,
            stringsAsFactors = FALSE
        )
    dismc <- dismc[order(dismc$class),]
    for (i in (1:n_cla)) {
        points(mu[i, 1], mu[i, 2])
    }
    points <- list()
    hulls <- list()
    for (i in (1:n_cla)) {
        points[[i]] <- cbind(dismc$x1[dismc$class == nam[i]],
                             dismc$x2[dismc$class == nam[i]])
        hulls[[i]] <- concaveman::concaveman(points[[i]], concavity = 10e20)
        lines(hulls[[i]])
    }
    return(list(pi, cov, mu))
}
lda2(X)
rnorm
