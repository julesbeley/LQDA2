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
                 x2 = rnorm(500, mean = 6, sd = 1),
                 class = "orange")
X <- rbind(X1, X2, X3)
plot(X$x1[X$class == "black"], X$x2[X$class == "black"], xlim = c(-9, 5), ylim = c(-1, 12))
points(X$x1[X$class == "white"], X$x2[X$class == "white"], col = "red")

lda2 <- function(value, data) {
    if (is.numeric(data[, 1] & is.numeric(data[, 2]))) {
        names(data) <- c("x1", "x2", "class")
    }
    tab <- table(data$class)
    n_cla <- dim(tab)
    nam <- names(tab)
    pi <- c()
    mu <- matrix(nrow = 2, ncol = n_cla)
    for (i in (1:n_cla)) {
        pi[i] <- tab[i] / dim(data)[1]
        mu[1, i] <- mean(data$x1[data$class == nam[i]]) # mean matrix: (class1:x1 ... classn:x1)
        mu[2, i] <- mean(data$x2[data$class == nam[i]]) #              (class1:x2 ... classn:x2)
    }
    val <- list()
    dev <- list()
    means <- matrix(nrow = n_cla, ncol = 2)
    mul <- list()
    ccov <- list()
    for (i in (1:n_cla)) {
        val[[i]] <- matrix(nrow = tab[i], ncol = 2)
        val[[i]] <- cbind(data$x1[data$class == nam[i]], data$x2[data$class == nam[i]])
        dev[[i]] <- matrix(nrow = tab[i], ncol = 2)
        means[i,] <- cbind(mean(data$x1[data$class == nam[i]]), mean(data$x2[data$class == nam[i]])) # redundant with mu!
        for (j in (1:tab[i])) {
            dev[[i]][j,] <- val[[i]][j, ] - means[i, ]
        }
        mul[[i]] <- matrix(nrow = 2, ncol = 2)
        ccov[[i]] <- matrix(nrow = 2, ncol = 2)
        mul[[i]] <- t(dev[[i]]) %*% dev[[i]] 
        ccov[[i]] <- mul[[i]] / (dim(data)[1] - n_cla) 
    }
    dis <- c()
    exp <- c()
    cov <- Reduce('+', ccov)
    inv <- solve(cov)
    for (i in (1:n_cla)) {
        dis[i] <- t(value) %*% inv %*% mu[, i] - (0.5 %*% t(mu[, i]) %*% inv %*% mu[, i]) + log(pi[i])
    }
    prob <- c()
    prob <- exp(dis) / sum(exp(dis))
    names(dis) <- nam
    class <- names(dis)[dis == max(dis)]
    score <- max(dis)
    proba <- max(prob)
    out <- as.data.frame(cbind(class, score, proba))
    return(out)
}
lda2(c(-2,5), X)
