# EMPTY HULLS BUG
# replace values near bounding box with actual bounding box values
# develop qda function

rm(list = ls())
Xlist <- list()
classnames <- c("white", "black", "blue", "red", "green", "orange", "purple", "brown", "car", "truck",
                "limo", "coke", "soda")
length <- runif(13, 100, 500)
for (i in (1:runif(1, 2, 13))) {
      Xlist[[i]] <- data.frame(
            x1 = rnorm(length[i], mean = runif(1, -10, 10), sd = runif(1, 0.5, 4)),
            x2 = rnorm(length[i], mean = runif(1, -10, 10), sd = runif(1, 0.5, 4)),
            class = classnames[i]
      )
}
do.call("rbind", Xlist) -> X

lda2 <- function(data, k = 4) {
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
      box <- data.frame(x1 = c(min(minx[, 1]),
                               rep.int(max(maxx[, 1]), 2),
                               min(minx[, 1])),
                        x2 = c(rep.int(min(minx[, 2]), 2),
                               rep.int(max(maxx[, 2]), 2)))
      if (all(pi %in% rep.int(pi[1], length(pi)))) {
            library(ggvoronoi)
            voronoi_polygon(data = as.data.frame(mu),
                            x = "x1", 
                            y = "x2",
                            outline = box) -> raw
            rm(minx, maxx)
            hulls <- list()
            for (i in (1:n_cla)) {
                  hulls[[i]] <- as.data.frame(raw@polygons[[i]]@Polygons[[1]]@coords)
                  names(hulls[[i]]) <- c("x1", "x2")
            }
      }
      else { 
            library(sp)
            runif <- cbind(runif(1000, min = min(minx[, 1]), max = max(maxx[, 1])), 
                           runif(1000, min = min(minx[, 2]), max = max(maxx[, 2])))
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
                  x1 = runif[, 1], 
                  x2 = runif[, 2],
                  class = classmc,
                  stringsAsFactors = FALSE
            )
            dismc <- dismc[order(dismc$class), ]
            test <- 0
            for (i in (1:n_cla)) {
                  if (length(dismc$x1[dismc$class %in% nam[i]]) == 0) test <- test + 1
            }
            if (test != 0) {
                  warning("Class missing: running patch")
                  add <- list()
                  for (i in (1:n_cla)) { 
                        if (length(dismc$x1[dismc$class %in% nam[i]]) == 0) {
                              cbind(rnorm(n = 300,
                                          mean = mu[i, 1],
                                          sd = sd(data$x1[data$class %in% nam[i]])),
                                    rnorm(n = 300,
                                          mean = mu[i, 2],
                                          sd = sd(data$x2[data$class %in% nam[i]])
                                    )) -> add[[i]]
                        }
                  }
                  do.call("rbind", add) -> add
                  dismc3 <- matrix(nrow = nrow(add), ncol = n_cla)
                  for (i in (1:n_cla)) {
                        for (h in (1:nrow(add))) {
                              dismc3[h, i] <- t(add[h, ]) %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
                        }
                  }
                  colnames(dismc3) <- nam
                  classmc3 <- c()
                  for (h in (1:300)) {
                        classmc3[h] <- names(dismc3[h, ])[dismc3[h, ] %in% max(dismc3[h, ])]
                  }
                  dismc3 <- data.frame(
                        x1 = add[, 1],
                        x2 = add[, 2],
                        class = classmc3
                  )
                  rbind(dismc, dismc3) -> dismc
                  dismc <- dismc[order(dismc$class), ]
            }
            mtosp <- function(m) SpatialPolygons(list(Polygons(list(Polygon(m)), 1)))
            box <- mtosp(box)
            hulls <- list()
            marker <- c()
            for (i in (1:n_cla)) { 
                  hulls[[i]] <- try(as.data.frame(
                        concaveman::concaveman(cbind(dismc$x1[dismc$class %in% nam[i]],
                                                     dismc$x2[dismc$class %in% nam[i]]),
                                               concavity = 20)),
                        silent = TRUE)
                  hulls[[i]] <- suppressWarnings(try(mtosp(hulls[[i]]), silent = TRUE)) 
                  if (class(hulls[[i]]) %in% "try-error") { # error message: class i too small
                        stop <- TRUE
                        marker[i] <- nam[i]
                  }
            }
            marker <- marker[!is.na(marker)]
            if (isTRUE(stop)) stop(paste("Class", marker, "is too small to be approximated"))
            for (j in (1:k)) {
                  for (i in (1:n_cla)) {
                        suppressWarnings(hulls[[i]] <- mtosp(hulls[[i]]))
                  }
                  rgeos::gDifference(box, hulls[[1]]) -> diff
                  for (i in (2:n_cla)) {
                        rgeos::gDifference(diff, hulls[[i]]) -> diff
                  }
                  spsample(diff, n = 2000, "random") -> sample
                  sample@coords -> sample
                  dismc2 <- matrix(nrow = 2000, ncol = n_cla)
                  for (h in (1:2000)) {
                        for (i in (1:n_cla)) {
                              dismc2[h, i] <- sample[h, ] %*% inv %*% mu[i, ] - 0.5 %*% t(mu[i, ]) %*% inv %*% mu[i, ] + log(pi[i])
                        }
                  }
                  colnames(dismc2) <- nam
                  classmc2 <- c()
                  for (h in (1:2000)) {
                        classmc2[h] <- names(dismc2[h, ])[dismc2[h, ] %in% max(dismc2[h, ])]
                  }
                  dismc2 <- data.frame(
                        x1 = sample[, 1], 
                        x2 = sample[, 2],
                        class = classmc2,
                        stringsAsFactors = FALSE
                  ) 
                  dismc2 <- dismc2[order(dismc2$class), ]
                  for (i in (1:n_cla)) {
                        hulls[[i]] <- as.data.frame(
                              concaveman::concaveman(cbind(dismc2$x1[dismc2$class %in% nam[i]],
                                                           dismc2$x2[dismc2$class %in% nam[i]]),
                                                     concavity = 20))
                        names(hulls[[i]]) <- c("x1", "x2")
                  }
                  rm(dismc2)
            }
      }
      list(pi, cov, mu, hulls) -> out
      names(out) <- c("Prior probabilities", 
                      "Covariance matrix", 
                      "Class means",
                      "Decision boundaries")
      return(out)
}
lda2(X) -> t

col <- heat.colors(length(table(X$class)))
ggplot() +
      geom_polygon(data = t$`Decision boundaries`[[1]],
                aes(x = x1, y = x2), fill = col[1], col = "black") -> r
for (i in (2:length(t$`Decision boundaries`))) {
      r + geom_polygon(data = t$`Decision boundaries`[[i]],
                    aes(x = x1, y = x2),
                    fill = col[i],
                    col = "black") + theme_void() -> r
}
r
table(X$class)
