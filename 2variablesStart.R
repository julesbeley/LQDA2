# generate random points on the fringes of the normal distribution + near bounding box
# develop qda function

rm(list = ls())
Xlist <- list()
classnames <- c("white", "black", "blue", "red", "green", "orange", "purple", "brown")
length <- runif(8, 100, 500)
for (i in (1:runif(1, 3, 8))) {
    Xlist[[i]] <- data.frame(
        x1 = rnorm(length[i], mean = runif(1, -10, 10), sd = 1),
        x2 = rnorm(length[i], mean = runif(1, -10, 10), sd = 1),
        class = classnames[i]
    )
}
do.call("rbind", Xlist) -> X

lda2 <- function(data, nmc = 200000, palette = heat.colors, ...) {
    library(ggplot2)
    library(dplyr)
    if (is.numeric(data[, 1] && is.numeric(data[, 2]))) {
        names(data) <- c("x1", "x2", "class")
    }
    tab <- table(data$class)
    n_cla <- dim(tab)
    nam <- names(tab)
    pi <- c()
    mu <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        pi[i] <- tab[i] / dim(data)[1]
        mu[i, 1] <- mean(data$x1[data$class %in% nam[i]]) 
        mu[i, 2] <- mean(data$x2[data$class %in% nam[i]])
    }
    names(pi) <- nam
    colnames(mu) <- c("x1", "x2")
    col <- palette(n = n_cla, alpha = 0.7)
    val <- list()
    dev <- list()
    mul <- list()
    ccov <- list()
    for (i in (1:n_cla)) {
        val[[i]] <- matrix(nrow = tab[i], ncol = 2)
        val[[i]] <- cbind(data$x1[data$class %in% nam[i]], data$x2[data$class %in% nam[i]])
        dev[[i]] <- matrix(nrow = tab[i], ncol = 2)
        for (j in (1:tab[i])) {
            dev[[i]][j, ] <- val[[i]][j, ] - mu[i, ]
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
    runifx1 <- runif(nmc,
                     min = min(minx[, 1]),
                     max = max(maxx[, 1]))
    runifx2 <- runif(nmc,
                     min = min(minx[, 2]),
                     max = max(maxx[, 2]))
    runif <- cbind(runifx1, runifx2)
    dismc <- matrix(nrow = nmc, ncol = n_cla)
    for (h in (1:nmc)) {
        for (i in (1:n_cla)) {
            dismc[h, i] <- t(runif[h, ]) %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
        }
    }
    colnames(dismc) <- nam
    rownames(mu) <- nam
    classmc <- c()
    for (h in (1:nmc)) {
        classmc[h] <- names(dismc[h, ])[dismc[h, ] %in% max(dismc[h, ])]
    }
    dismc <- data.frame(
        x1 = runifx1, 
        x2 = runifx2,
        class = classmc,
        stringsAsFactors = FALSE
        )
    dismc <- dismc[order(dismc$class), ]
    points <- list()
    hulls <- list()
    ggplot() + 
        geom_point(data = subset(data, class %in% nam[1]),
                   aes(x = x1, y = x2),
                   col = col[1]) +
        scale_x_continuous(limits = c(1.1 * min(minx[, 1]), 1.1 * max(maxx[, 1]))) +
        scale_y_continuous(limits = c(1.1 * min(minx[, 2]), 1.1 * max(maxx[, 2]))) +
        xlab("X1") +
        ylab("X2") +
        labs(title = "Observations and predicted decision boundaries") -> g
    for (i in (2:n_cla)) {
        g + geom_point(data = subset(data, class %in% nam[i]),
                       aes(x = x1, y = x2),
                       col = col[i],
                       alpha = 0.4) + 
            theme(plot.title = element_text(hjust = 0.5),
                  panel.background = element_blank(),
                  ...) -> g
    }
    if (all(pi %in% rep.int(pi[1], length(pi)))) {
        library(ggvoronoi)
        box <- data.frame(x1 = c(min(minx[, 1]),
                                 rep.int(max(maxx[, 1]), 2),
                                 min(minx[, 1])),
                          x2 = c(rep.int(min(minx[, 2]), 2),
                                 rep.int(max(maxx[, 2]), 2)))
        g + stat_voronoi(data = as.data.frame(mu),
                         aes(x = x1, y = x2),
                         geom = "path",
                         outline = box,
                         size = 1,
                         col = "grey30") -> g
    }
    else { 
        for (i in (1:n_cla)) {
        points[[i]] <- cbind(dismc$x1[dismc$class %in% nam[i]],
                             dismc$x2[dismc$class %in% nam[i]])
        hulls[[i]] <- concaveman::concaveman(points[[i]], concavity = 10e20)
        as.data.frame(hulls[[i]]) -> hulls[[i]]
        g + geom_path(data = hulls[[i]],
                      aes(x = V1, y = V2),
                      size = 1,
                      col = "grey30") -> g 
        }
    }
    print(g)
    list(pi, cov, mu) -> out
    names(out) <- c("Prior probabilities", 
                    "Covariance matrix", 
                    "Class means")
    return(out)
}

lda2(data = X, nmc = 200000, palette = heat.colors)

