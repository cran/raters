concordance <- function(db, test = "Default", B = 1000, alpha = 0.05) {
  db <- as.matrix(db)
  db2 <- db * (db - 1)
  R <- nrow(db)
  C <- ncol(db)
  sumcol <- colSums(db)
  sumrow <- rowSums(db)
  tot <- sum(db)
  vet <- list()
  pij <- matrix(NA, nrow = R, ncol = C)

  pij <- db2 / (sumrow * (sumrow - 1))

  pi <- rowSums(pij)
  p <- mean(pi)
  pj <- sumcol / tot
  pj2 <- pj^2
  pe <- sum(pj2)
  fleiss.kappa <- (p - pe) / (1 - pe)
  s <- (C * p - 1) / (C - 1)
  s.boot <- c()
  pi.v.boot <- replicate(B, pi.boot <- sample(pi,
                                              size = R,
                                              replace = TRUE
  ))
  p.boot <- apply(pi.v.boot, 2, mean)
  s.boot <- (C * p.boot - 1) / (C - 1)
  s.boot.ci <- quantile(s.boot, probs = c(alpha / 2, 1 - alpha / 2))

  Default <- function(vet) {
    s.vet <- c(s, s.boot.ci)
    names(s.vet) <- c("S", "LCL", "UCL")
    vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
    cat(paste("Inter-rater Agreement"), "\n")
    message("Fleiss kappa confidence intervals and test not estimated since the number of raters is not the same among evaluated subjects (n)")
    vet
  }

  Normal <- function(db) {
    stat.test <- s * R * sqrt((C - 1) / (2 * sum(1 / (sumrow * (sumrow - 1)))))
    bin <- stat.test > qnorm(1 - alpha)
    pvalue <- 1 - pnorm(stat.test)
    s.vet <- c(s, s.boot.ci, pvalue)
    names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
    vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
    cat(paste("Inter-rater Agreement"), "\n")
    message("Fleiss kappa confidence intervals and test not estimated since the number of raters is not the same among evaluated subjects (n)")
    vet
  }

  Chisq <- function(db)
    {stop("Chi-squared test approximation can be used only when total number of raters is the same among evaluated subjects (n)")}

  MC <- function(db) {

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

    pij.mc <- lapply(matrix.mc,
                            function(list) {
                              list * (list - 1) / (sumrow * (sumrow - 1))
                            })

    pi.mc <- lapply(pij.mc, function(list) {apply(list, 1, sum)})

    p.mc <- do.call(rbind, lapply(pi.mc,mean))
    s.mc <- ((p.mc * C) - 1) / (C - 1)
    crit.s.mc <- quantile(s.mc, 0.95)
    binary <- (s.mc >= s)
    pvalue <- sum(binary) / B
    s.vet <- c(s, s.boot.ci, pvalue)
    names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
    vet <- list(Fleiss = fleiss.kappa, Statistic = s.vet)
    cat(paste("Inter-rater Agreement"), "\n")
    message("Fleiss kappa confidence intervals and test not estimated since the number of raters is not the same among evaluated subjects (n)")
    vet
  }

  if (length(unique(rowSums(db))) == 1) {
    std.kappaj <- sqrt(2 / (R * (sumrow[1] * (sumrow[1] - 1))))
    pj.1 <- 1 - pj
    pj.2 <- 1 - 2 * pj
    std.num <- sqrt((sum(pj * pj.1))^2 - sum(pj * pj.1 * pj.2))
    std.den <- sum(pj * pj.1)
    std.adj <- std.num / std.den
    std.kappa <- std.kappaj * std.adj
    ci.low <- fleiss.kappa - qnorm(1 - alpha / 2) * std.kappa
    ci.upp <- fleiss.kappa + qnorm(1 - alpha / 2) * std.kappa
    z.fleiss <- fleiss.kappa / std.kappa
    pvalue.fleiss <- 1 - pnorm(z.fleiss)
    if (ci.upp > 1 | ci.low < (-1 / (sumrow[1] - 1))) {
      ci.low <- (1 / (sumrow[1] - 1))
      ci.upp <- 1
    }

    Default2 <- function(vet) {
      kappa.vet <- c(
        fleiss.kappa,
        ci.low,
        ci.upp,
        std.kappa
      )
      names(kappa.vet) <- c(
        "Kappa",
        "LCL",
        "UCL",
        "Std.Error"
      )

      s.vet <- c(s, s.boot.ci)
      names(s.vet) <- c("S", "LCL", "UCL")
      vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
      cat(paste("Inter-rater Agreement"), "\n")
      vet
    }

    Chisq2 <- function(db) {
      stat.test <- R * (C - 1) * ((sumrow -1)[1] * s + 1)
      bin <- stat.test > qchisq(alpha, df = R *(C - 1))
      pvalue <- 1 - pchisq(stat.test, df = R *(C - 1))
      kappa.vet <- c(
        fleiss.kappa,
        ci.low,
        ci.upp,
        std.kappa,
        z.fleiss,
        pvalue.fleiss
      )
      names(kappa.vet) <- c(
        "Kappa",
        "LCL",
        "UCL",
        "Std.Error",
        "Z value",
        "Pr(>|z|)"
      )
      s.vet <- c(s, s.boot.ci, pvalue)
      names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
      vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
      cat(paste("Inter-rater Agreement"), "\n")
      vet
    }

    Normal2 <- function(db) {
      stat.test <- s * R * sqrt((C - 1) / (2 * sum(1 / (sumrow * (sumrow - 1)))))
      bin <- stat.test > qnorm(1 - alpha)
      pvalue <- 1 - pnorm(stat.test)
      kappa.vet <- c(
        fleiss.kappa,
        ci.low,
        ci.upp,
        std.kappa,
        z.fleiss,
        pvalue.fleiss
      )
      names(kappa.vet) <- c(
        "Kappa",
        "LCL",
        "UCL",
        "Std.Error",
        "Z value",
        "Pr(>|z|)"
      )
      s.vet <- c(s, s.boot.ci, pvalue)
      names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
      vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
      cat(paste("Inter-rater Agreement"), "\n")
      vet
    }

    MC2 <- function(db) {

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

      pij.mc <- lapply(matrix.mc,
                       function(list) {
                         list * (list - 1) / (sumrow * (sumrow - 1))
                       })

      pi.mc <- lapply(pij.mc, function(list) {apply(list, 1, sum)})

      p.mc <- do.call(rbind, lapply(pi.mc,mean))
      s.mc <- ((p.mc * C) - 1) / (C - 1)
      crit.s.mc <- quantile(s.mc, 0.95)
      binary <- (s.mc >= s)
      pvalue <- sum(binary) / B

      kappa.vet <- c(
        fleiss.kappa,
        ci.low,
        ci.upp,
        std.kappa,
        z.fleiss,
        pvalue.fleiss
      )

      names(kappa.vet) <- c(
        "Kappa",
        "LCL",
        "UCL",
        "Std.Error",
        "Z value",
        "Pr(>|z|)"
      )
      s.vet <- c(s, s.boot.ci, pvalue)
      names(s.vet) <- c("S", "LCL", "UCL", "pvalue")
      vet <- list(Fleiss = kappa.vet, Statistic = s.vet)
      cat(paste("Inter-rater Agreement"), "\n")
      vet
    }

    switch(test,
           Normal = Normal2(db),
           MC = MC2(db),
           Chisq = Chisq2(db),
           Default = Default2(vet)
    )
  }
  else {
    switch(test,
           Normal = Normal(db),
           MC = MC(db),
           Chisq = Chisq(db),
           Default = Default(vet)
    )
  }
}




