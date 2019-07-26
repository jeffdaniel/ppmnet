test_that("(composite) likelihood calculations are correct", {

  # Pseudolikelihood
  Qp <- quadscheme(Xp)
  Qs <- quadscheme(Xs)
  Qg <- quadscheme(Xg)
  fit1p <- ppm(Qp ~ ., data = exdata)
  fit2p <- ppmnet(Qp, data = exdata, lambda.min.ratio = 1e-9)
  fit1s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5))
  fit2s <- ppmnet(Qs, data = exdata, interaction = Strauss(5),
                  lambda.min.ratio = 1e-9)
  fit1g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1))
  fit2g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1),
                  lambda.min.ratio = 1e-9)
  expect_equal(as.numeric(logLik(fit1p)),
               as.numeric(logLik(fit2p)[length(fit2p$lambda)]),
               tolerance = 0.001)
  expect_equal(as.numeric(logLik(fit1s)),
               as.numeric(logLik(fit2s)[length(fit2s$lambda)]),
               tolerance = 0.001)
  expect_equal(as.numeric(logLik(fit1g)),
               as.numeric(logLik(fit2g)[length(fit2g$lambda)]),
               tolerance = 0.001)

  # Logistic Composite Likelihood
  Qp <- quadscheme.logi(Xp)
  Qs <- quadscheme.logi(Xs)
  Qg <- quadscheme.logi(Xg)
  fit1p <- ppm(Qp ~ ., data = exdata, method = "logi")
  fit2p <- ppmnet(Qp, data = exdata, lambda.min.ratio = 1e-9, method = "logi")
  fit1s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5),
               method = "logi")
  fit2s <- ppmnet(Qs, data = exdata, interaction = Strauss(5),
                  method = "logi", lambda.min.ratio = 1e-9)
  fit1g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1),
               method = "logi")
  fit2g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1),
                  method = "logi", lambda.min.ratio = 1e-9)
  expect_equal(as.numeric(logLik(fit1p)),
               as.numeric(logLik(fit2p)[length(fit2p$lambda)]),
               tolerance = 0.001)
  expect_equal(as.numeric(logLik(fit1s)),
               as.numeric(logLik(fit2s)[length(fit2s$lambda)]),
               tolerance = 0.001)
  expect_equal(as.numeric(logLik(fit1g)),
               as.numeric(logLik(fit2g)[length(fit2g$lambda)]),
               tolerance = 0.001)

})

test_that("effective degrees of freedom are correct", {

  # Pseudolikelihood
  Qs <- quadscheme(Xs)
  fit1s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5))
  fit2s <- ppmnet(Qs, data = exdata, interaction = Strauss(5),
                  nlambda = 10, lambda.min.ratio = 1e-9)
  V <- vcov(fit1s, what = "internals")
  J <- V$Sigma
  H <- V$A1
  JiH <- solve(H, J)
  expect_equal(sum(diag(JiH)), select(fit2s)$edf[length(fit2s$lambda)],
               tolerance = 1e-3)

  Qg <- quadscheme(Xg)
  fit1g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1))
  fit2g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1),
                  nlambda = 10, lambda.min.ratio = 1e-9)
  V <- vcov(fit1g, what = "internals")
  J <- V$Sigma
  H <- V$A1
  JiH <- solve(H, J)
  expect_equal(sum(diag(JiH)), select(fit2g)$edf[length(fit2g$lambda)],
               tolerance = 1e-3)

  # Logistic Composite Likelihood
  Qs <- quadscheme.logi(Xs)
  fit1s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5),
               method = "logi")
  fit2s <- ppmnet(Qs, data = exdata, interaction = Strauss(5), method = "logi",
                  nlambda = 10, lambda.min.ratio = 1e-9)
  V <- vcov(fit1s, what = "internals")
  J <- V$Sigma1log + V$Sigma2log
  H <- V$Slog
  JiH <- solve(H, J)
  expect_equal(sum(diag(JiH)), select(fit2s)$edf[length(fit2s$lambda)],
               tolerance = 0.1)

  Qg <- quadscheme.logi(Xg)
  fit1g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1),
               method = "logi")
  fit2g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1), method = "logi",
                  nlambda = 10, lambda.min.ratio = 1e-9)
  V <- vcov(fit1g, what = "internals")
  J <- V$Sigma1log + V$Sigma2log
  H <- V$Slog
  JiH <- solve(H, J)
  expect_equal(sum(diag(JiH)), select(fit2g)$edf[length(fit2g$lambda)],
               tolerance = 0.1)


})
