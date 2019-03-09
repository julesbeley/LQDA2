# GENERATE RANDOM POINTS FURTHER THAN BOUNDARY WITH FEWER RANDOM POINTS!!!!!!!
# separate graphical part from computation of decision boundary path
# develop qda function

rm(list = ls())
Xlist <- list()
classnames <- c("white", "black", "blue", "red", "green", "orange", "purple", "brown")
length <- runif(8, 100, 500)
for (i in (1:runif(1, 3, 8))) {
    Xlist[[i]] <- data.frame(
        x1 = rnorm(length[i], mean = runif(1, -10, 10), sd = runif(1, 0.5, 4)),
        x2 = rnorm(length[i], mean = runif(1, -10, 10), sd = runif(1, 0.5, 4)),
        class = classnames[i]
    )
}
do.call("rbind", Xlist) -> X

lda2 <- function(data, palette = heat.colors, ...) {
    library(ggplot2)
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
        val[[i]] <- cbind(data$x1[data$class %in% nam[i]], 
                          data$x2[data$class %in% nam[i]])
        dev[[i]] <- matrix(nrow = tab[i], ncol = 2)
        for (j in (1:tab[i])) {
            dev[[i]][j, ] <- val[[i]][j, ] - mu[i, ]
        }
        ccov[[i]] <- matrix(nrow = 2, ncol = 2)
        ccov[[i]] <- t(dev[[i]]) %*% dev[[i]] / (dim(data)[1] - n_cla)
    }
    dis <- c()
    exp <- c()
    cov <- Reduce('+', ccov)
    rm(val, dev, ccov)
    rownames(cov) <- c("x1", "x2")
    colnames(cov) <- c("x1", "x2")
    inv <- solve(cov)
    maxx <- matrix(nrow = n_cla, ncol = 2)
    minx <- matrix(nrow = n_cla, ncol = 2)
    for (i in (1:n_cla)) {
        for (j in (1:2)) {
            maxx[i, j] <- mu[i, j] + 6 * sqrt(cov[j, j])
            minx[i, j] <- mu[i, j] - 6 * sqrt(cov[j, j])
        }
    }
    ggplot() + 
        geom_point(data = subset(data, class %in% nam[1]),
                   aes(x = x1, y = x2),
                   col = col[1]) +
        scale_x_continuous(limits = c(1.15 * min(minx[, 1]), 1.15 * max(maxx[, 1]))) +
        scale_y_continuous(limits = c(1.15 * min(minx[, 2]), 1.15 * max(maxx[, 2]))) +
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
    box <- data.frame(x1 = c(min(minx[, 1]),
                             rep.int(max(maxx[, 1]), 2),
                             min(minx[, 1])),
                      x2 = c(rep.int(min(minx[, 2]), 2),
                             rep.int(max(maxx[, 2]), 2)))
    if (all(pi %in% rep.int(pi[1], length(pi)))) {
        g + ggvoronoi::stat_voronoi(data = as.data.frame(mu),
                                    aes(x = x1, y = x2),
                                    geom = "path",
                                    outline = box,
                                    size = 1,
                                    col = "grey30") -> g
        rm(minx, maxx)
    }
    else { 
        library(sp)
        runifx1 <- runif(1000,
                         min = min(minx[, 1]),
                         max = max(maxx[, 1]))
        runifx2 <- runif(1000,
                         min = min(minx[, 2]),
                         max = max(maxx[, 2]))
        runif <- cbind(runifx1, runifx2)
        dismc <- matrix(nrow = 1000, ncol = n_cla)
        for (h in (1:1000)) {
            for (i in (1:n_cla)) {
                dismc[h, i] <- t(runif[h, ]) %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
            }
        }
        colnames(dismc) <- nam
        rownames(mu) <- nam
        classmc <- c()
        for (h in (1:1000)) {
            classmc[h] <- names(dismc[h, ])[dismc[h, ] %in% max(dismc[h, ])]
        }
        dismc <- data.frame(
            x1 = runifx1, 
            x2 = runifx2,
            class = classmc,
            stringsAsFactors = FALSE
        )
        rm(runifx1, runifx2)
        dismc <- dismc[order(dismc$class), ]
        hulls <- list()
        mtosp <- function(m) SpatialPolygons(list(Polygons(list(Polygon(m)), 1)))
        box <- mtosp(box)
        for (i in (1:n_cla)) {
            hulls[[i]] <- as.data.frame(
                concaveman::concaveman(cbind(dismc$x1[dismc$class %in% nam[i]],
                                             dismc$x2[dismc$class %in% nam[i]]),
                                       concavity = 50))
            hulls[[i]] <- mtosp(hulls[[i]])
        }
        rgeos::gDifference(box, hulls[[1]]) -> diff
        for (i in (2:n_cla)) {
            rgeos::gDifference(diff, hulls[[i]]) -> diff
        }
    }
    print(g)
    list(pi, cov, mu, diff) -> out
    names(out) <- c("Prior probabilities", 
                    "Covariance matrix", 
                    "Class means",
                    "Hulls")
    return(out)
}
lda2(X) -> r

plot(r$Hulls)
spsample(r$Hulls, n = 200, "random")
points(spsample(r$Hulls, n = 200, "random"))


# plotting section

for (i in (1:n_cla)) {
    g + geom_path(data = hulls[[i]],
                  aes(x = V1, y = V2),
                  size = 1,
                  col = "grey30") -> g 
}

# sampling in area outside of a polygon!

library(sp)
library(rgeos)
cbind(c(0, 3, 3, 0, 0), c(0, 0, 3, 3, 0)) -> pol
cbind(c(0, 5, 5, 0, 0), c(0, 0, 5, 5, 0)) -> pol2
cbind(c(1, 2, 2, 1, 1), c(1, 1, 2, 2, 1)) -> pol3
sps <- SpatialPolygons(list(Polygons(list(Polygon(pol)), 1)))
sps2 <- SpatialPolygons(list(Polygons(list(Polygon(pol2)), 1)))
sps3 <- SpatialPolygons(list(Polygons(list(Polygon(pol3)), 1)))

gDifference(sps2, sps3) -> diff
plot(diff)
points(spsample(diff, n = 2000, "random"))


# convert back to data.frame for testing
sps@polygons[[1]]@Polygons[[1]]@coords -> d
diff@polygons[[1]]@Polygons[[1]]@coords 


# gif test

# for (i in seq(from = 50, to = 10050, by = 100)) {
#    png(print(paste("./plot", i , ".png", sep = "")))
#    lda2(X, nmc = i)
#    dev.off()
#}