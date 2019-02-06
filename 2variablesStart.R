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
X <- rbind(X0, X1, X2, X3) 

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
    slopes <- matrix(nrow = n_cla, ncol = n_cla) # slope of mean lines between class i and j
    intercepts <- matrix(nrow = n_cla, ncol = n_cla) # intercept of mean lines "
    for (i in (1:n_cla)) {
        for (j in (i:n_cla)) {
            slopes[i,j] <- (mu[j,2] - mu[i,2]) / (mu[j,1] - mu[i,1])
            intercepts[i,j] <- mu[i,2] - slopes[i,j] * mu[i,1]
        }
    }
    x <- seq(mu[1,1], mu[2,1], length = 100)
    y <- slopes[1,2] * x + intercepts[1,2]
    xy <- cbind(x,y)
    dis <- matrix(nrow = 100, ncol = n_cla)
    for (j in (1:100)) {
        for (i in (1:n_cla))
            dis[j,i] <- t(xy[j,]) %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
    }
    colnames(dis) <- nam
    classtest <- c()
    for (j in (1:100)) {
        classtest[j] <- colnames(dis)[dis[j,] == max(dis[j,])]
    }
    dis <- as.data.frame(cbind(dis,classtest)) # scores and class for every point on the means line
    if (length(table(dis$classtest)) == 2) {
        for (j in (1:99)) {
            if (dis$classtest[j] != dis$classtest[j+1]) {
                threshold <- j
            }
        }
    }
    bounds1 <- data.frame(x1 = xy[threshold,1], x2 = xy[threshold,2])
    bounds2 <- data.frame(x1 = xy[threshold + 1,1], x2 = xy[threshold + 1,2])
    bounds <- rbind(bounds1, bounds2)
    rownames(bounds) <- NULL
    return(list(pi, cov, mu, bounds))
}
lda2(X)
