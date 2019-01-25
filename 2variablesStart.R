X1 <-
    data.frame(
        x1 = rnorm(1000, mean = 0, sd = 1),
        x2 = rnorm(1000, mean = 3, sd = 1),
        class = "black"
    )
X2 <-
    data.frame(
        x1 = rnorm(2000, mean = -4, sd = 1),
        x2 = rnorm(2000, mean = 8, sd = 1),
        class = "white"
    )
X <- rbind(X1, X2)
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
        mu[1, i] <- mean(data$x1[data$class == nam[i]]) # mean matrix: (class1:x1 ... class2:x1)
        mu[2, i] <- mean(data$x2[data$class == nam[i]]) #              (class1:x2 ... class2:x2)
    }
    values <- list()
    deviation <- list()
    means <- matrix(nrow = n_cla, ncol = 2)
    multiplied <- list()
    ccovariance <- list()
    for (i in (1:n_cla)) {
        values[[i]] <- matrix(nrow = tab[i], ncol = 2)
        values[[i]] <- cbind(data$x1[data$class == nam[i]], data$x2[data$class == nam[i]])
        deviation[[i]] <- matrix(nrow = tab[i], ncol = 2)  
        means[i,] <- cbind(mean(data$x1[data$class == nam[i]]), mean(data$x2[data$class == nam[i]]))
        for (j in (1:tab[i])) {
            deviation[[i]][j,] <- values[[i]][j, ] - means[i, ]
        }
        multiplied[[i]] <- matrix(nrow = 2, ncol = 2)
        ccovariance[[i]] <- matrix(nrow = 2, ncol = 2)
        multiplied[[i]] <- t(deviation[[i]]) %*% deviation[[i]]
        ccovariance[[i]] <- multiplied[[i]] / (dim(data)[1] - n_cla)  
    }
    covariance <- Reduce('+', ccovariance)
    return(covariance)
}
lda2(1, X)

