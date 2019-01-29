X0 <-
    data.frame(
        x1 = rnorm(400, mean = -12, sd = 1),
        x2 = rnorm(400, mean = 7, sd = 1),
        class = "black"
    )
X1 <-
    data.frame(
        x1 = rnorm(200, mean = -7, sd = 1),
        x2 = rnorm(200, mean = 3, sd = 1),
        class = "black"
    )
X2 <-
    data.frame(
        x1 = rnorm(200, mean = -4, sd = 1),
        x2 = rnorm(200, mean = 6, sd = 1),
        class = "white"
    )
X3 <- data.frame(
    x1 = rnorm(200, mean = -1, sd = 1),
    x2 = rnorm(200, mean = -5, sd = 1),
    class = "orange"
)
X4 <- data.frame(
    x1 = rnorm(300, mean = 2, sd = 1),
    x2 = rnorm(300, mean = -2, sd = 1),
    class = "red"
)
X5 <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 1),
    x2 = rnorm(200, mean = 15, sd = 1),
    class = "blue"
)
X <- rbind(X0,X1, X2, X3, X4, X5)

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
            dev[[i]][j,] <- val[[i]][j, ] - mu[i, ]
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
    col <- heat.colors(n = n_cla)
    runifx1 <- runif(100000, min = min(minx[, 1]), max = max(maxx[, 1]))
    runifx2 <- runif(100000, min = min(minx[, 2]), max = max(maxx[, 2]))
    runif <- cbind(runifx1, runifx2)
    dismc <- matrix(nrow = 100000, ncol = n_cla)
    for (h in (1:100000)) {
        for (i in (1:n_cla)) {
            dismc[h, i] <-
                t(runif[h, ]) %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
        }
    }
    colnames(dismc) <- nam
    rownames(mu) <- nam
    classmc <- c()
    for (h in (1:100000)) {
        classmc[h] <- names(dismc[h, ])[dismc[h, ] == max(dismc[h, ])]
    }
    dismc <-
        data.frame(
            x1 = runifx1,
            x2 = runifx2,
            class = classmc,
            stringsAsFactors = FALSE
        )
    dismc <- dismc[order(dismc$class), ]
    plot(
        dismc$x1[dismc$class == nam[1]],
        dismc$x2[dismc$class == nam[1]],
        xlim = c(min(minx[, 1]), max(maxx[, 1])), # x1
        ylim = c(min(minx[, 2]), max(maxx[, 2])), # x2
        col = col[1],
        xlab = "X1",
        ylab = "X2",
        pch = 20,
        cex = 0.5
    )
    for (i in (2:n_cla)) {
        points(dismc$x1[dismc$class == nam[i]],
               dismc$x2[dismc$class == nam[i]], 
               col = col[i],
               pch = 20,
               cex = 0.5)
    }
    for (i in (1:n_cla)) {
        points(mu[i,1], mu[i,2],
               pch = )
    }
    for (i in (1:n_cla)) {
        for (j in (i:n_cla)) {
            segments(mu[i,1], mu[i,2], mu[j,1], mu[j,2])
        }
    }
    return(list(pi, cov, mu))
}
lda2(X)
# y = ax + b 
# a = (yB - yA) / (xB - xA)
# A c y so yA = aXA + b so
# b = yA - aXA

# if mean line doesn't cross right decision boundary, rotate slightly, if crosses two, chuck it

# OR (simpler algorithm) :
# branch out with two lines from the origin at 10 degrees (less?) on each side of average line and at 3/4 distance and
# retreat until function classifies to the original class (you have two points on the boundary so you have a boundary)
# THEN glide along this boundary until find new class - how far? how can you bound the search (to one side?)
