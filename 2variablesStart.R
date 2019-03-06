# convert to ggplot2
# create classifying function inside lda function 
# generate random points on the fringes of the normal distribution
# develop qda function

X1 <- data.frame(
    x1 = rnorm(100, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(100, mean = runif(1, -10, 10), sd = 1),
    class = "black"
    )
X2 <- data.frame(
    x1 = rnorm(200, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(200, mean = runif(1, -10, 10), sd = 1),
    class = "blue"
    )
X3 <- data.frame(
    x1 = rnorm(2000, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(2000, mean = runif(1, -10, 10), sd = 1),
    class = "orange"
)
X4 <- data.frame(
    x1 = rnorm(500, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(500, mean = runif(1, -10, 10), sd = 1),
    class = "red"
)
X5 <- data.frame(
    x1 = rnorm(600, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(600, mean = runif(1, -10, 10), sd = 1),
    class = "white"
)
X6 <- data.frame(
    x1 = rnorm(800, mean = runif(1, -10, 10), sd = 1),
    x2 = rnorm(800, mean = runif(1, -10, 10), sd = 1),
    class = "green"
)
X <- rbind(X1, X2, X3, X4, X5, X6)

lda2 <- function(data) {
    library(ggplot2)
    library(dplyr)
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
    maxx <- matrix(nrow = n_cla, ncol = 2)
    minx <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        for (j in (1:2)) {
            maxx[i, j] <- mu[i, j] + 4 * sqrt(cov[j, j])
            minx[i, j] <- mu[i, j] - 4 * sqrt(cov[j, j])
        }
    }
    col <- heat.colors(n = n_cla, alpha = 0.7)
    ggplot() + 
        geom_point(data = subset(data, class %in% nam[1]),
                   aes(x = x1, y = x2),
                   col = col[1]) +
        scale_x_continuous(limits = c(min(minx[, 1]), max(maxx[, 1]))) +
        scale_y_continuous(limits = c(min(minx[, 2]), max(maxx[, 2]))) +
        xlab("X1") +
        ylab("X2") +
        labs(title = "Observations and predicted decision boundaries") -> g
    for (i in (2:n_cla)) {
        g + geom_point(data = subset(data, class %in% nam[i]),
                       aes(x = x1, y = x2),
                       col = col[i]) + 
            theme(plot.title = element_text(hjust = 0.5)) -> g
    }
    runifx1 <- runif(200000,
                     min = 1.2 * min(minx[, 1]),
                     max = 1.2 * max(maxx[, 1]))
    runifx2 <- runif(200000,
                     min = 1.3 * min(minx[, 2]),
                     max = 1.3 * max(maxx[, 2]))
    runif <- cbind(runifx1, runifx2)
    dismc <- matrix(nrow = 200000, ncol = n_cla)
    for (h in (1:200000)) {
        for (i in (1:n_cla)) {
            dismc[h, i] <- t(runif[h,]) %*% inv %*% mu[i,] - 0.5 %*% t(mu[i,]) %*% inv %*% mu[i,] + log(pi[i])
        }
    }
    colnames(dismc) <- nam
    rownames(mu) <- nam
    classmc <- c()
    for (h in (1:200000)) {
        classmc[h] <- names(dismc[h,])[dismc[h,] == max(dismc[h,])]
    }
    dismc <- data.frame(
        x1 = runifx1, 
        x2 = runifx2,
        class = classmc,
        stringsAsFactors = FALSE
        )
    dismc <- dismc[order(dismc$class),]
    points <- list()
    hulls <- list()
    for (i in (1:n_cla)) {
        points[[i]] <- cbind(dismc$x1[dismc$class == nam[i]],
                             dismc$x2[dismc$class == nam[i]])
        hulls[[i]] <- concaveman::concaveman(points[[i]], concavity = 5)
        as.data.frame(hulls[[i]]) -> hulls[[i]]
    }
    for (i in (1:n_cla)) {
        for (j in (seq(1, dim(hulls[[i]])[1] - 1))) {
            if (abs(hulls[[i]][j, 1] - hulls[[i]][j + 1, 1]) < 0.02) {
                hulls[[i]][j,] <- c(NA, NA)
            }}
        g + geom_point(data = hulls[[i]], aes(x = V1, y = V2)) +
            geom_path(data = hulls[[i]], aes(x = V1, y = V2)) -> g
    }
    print(g)
    return(list(pi, cov, mu, hulls[[1]][1,1] - hulls[[1]][2,1] < 0.05))
}
lda2(X)  


# why the hell don't the hulls go all the way around????

# to test if priors are equal before applying ggvoronoi
x <- rep.int(10, 10)
x2 <- c(1, 2, 1, 1, 1, 10, 10)
x2 == rep.int(x2[1], length(x2))
test <- function(vector) {
    if (all(vector == rep.int(vector[1], length(vector)))) {
        print("same")
    }
    else {
        print("different")
    }
}
test(x2)

# "inverse" of a normal distribution to generate points far from centroid
