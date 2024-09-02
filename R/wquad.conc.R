wquad.conc <-
  function(db, test = "Default", B = 1000, alpha = 0.05) {
    db <- as.matrix(db)
    C <- ncol(db)
    R <- nrow(db)
    row.sum <- rowSums(db)
    db2 <- db^2
    row.sum2 <- rowSums(db2)

    w <- matrix(NA, nrow = C, ncol = C)
    for (j in 1:C) {
      for (k in 1:C) {
        w[j, k] <- 1 - ((abs(j - k))^2 / (C - 1)^2)
      }
    }

    ww <- rep(0, R)

    for (i in 1:R) {
      for (j in 1:(C - 1)) {
        l <- j + 1
        for (k in l:C) {
          ww[i] <- ww[i] + (db[i, j] * db[i, k] * w[j, k])
        }
      }
    }

    xi <- (0.5 * row.sum2) + ww - (0.5 * row.sum)
    pi <- (2 * xi) / (row.sum * (row.sum - 1))

    p.avg <- mean(pi)
    pe <- (1 / (C^2)) * sum(w)

    s.star <- (p.avg - pe) / (1 - pe)
    s.star

    pi.boot.w <- list()
    for (i in 1:B) {
      pi.boot.w[[i]] <- sample(pi, size = nrow(db), replace = TRUE)
    }
    p.boot.w <- do.call(rbind, lapply(pi.boot.w, mean))
    s.boot.w <- (p.boot.w - pe) / (1 - pe)
    s.boot.ci.w <- quantile(s.boot.w, probs = c(alpha / 2, 1 - alpha / 2))

    Default <- function(db) {
      s.res <- c(s.star, s.boot.ci.w)
      names(s.res) <- c("S*", "LCL", "UCL")
      s.res
    }

    MC <- function(db) {
      R <- nrow(db)
      C <- ncol(db)
      w.sum <- lapply(1:B, function(x) {rep(0, R)})
      xi.mc <- do.call(rbind, w.sum)
      pi.mc <- xi.mc

      matrix.mc <- lapply(1:B,
                          matrix,
                          data = NA,
                          nrow = R,
                          ncol = C)

      matrix.mc <- lapply(matrix.mc,
                          function(list) {
                            for (j in 1:R) {
                              list[j, ] <- t(rmultinom(1, size = rowSums(db)[j], prob = rep(1 / C, C)))
                            }
                            return(list)
                          })

      matrix.mc2 <- lapply(matrix.mc, function(x) {x^2})
      rs.mc <- do.call(rbind, lapply(matrix.mc, rowSums))
      rs2.mc <- do.call(rbind, lapply(matrix.mc2, rowSums))


      for (g in 1:B) {
        for (i in 1:R) {
          for (j in 1:(C - 1)) {
            l <- j + 1
            for (k in l:C) {
              w.sum[[g]][i] <- w.sum[[g]][i] + (matrix.mc[[g]][i, j] * matrix.mc[[g]][i, k] * w[j, k])
            }
          }
        }
      }

      w.sum <- do.call(rbind, w.sum)

      xi.mc <- (0.5 * rs2.mc)  + w.sum - (0.5 * rs.mc)
      pi.mc <- (2 * xi.mc) / (rs.mc *(rs.mc - 1))
      p.avg.mc <- rowMeans(pi.mc)
      pe.mc <- (1 / (C^2)) * sum(w)
      s.star.mc <- (p.avg.mc - pe.mc) / (1 - pe.mc)

      crit.s.mc <- quantile(s.star.mc, 0.95)
      binary <- s.star.mc >= s.star
      pvalue <- sum(binary) / B
      s.res <- c(s.star, s.boot.ci.w, pvalue)
      names(s.res) <- c("S*", "LCL", "UCL", "pvalue")
      s.res
    }

    switch(test,
           Default = Default(db),
           MC = MC(db)
    )
  }
