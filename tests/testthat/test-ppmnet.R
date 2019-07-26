test_that("ppmnet internals are consistent with ppm -- mpl", {

  # Poisson model
  Q <- quadscheme(bei)
  fit1 <- ppm(Q ~ ., data = bei.extra)
  fit2 <- ppmnet(Q, data = bei.extra)
  expect_identical(as.matrix(fit1$internal$glmdata[,-c(1,2,5)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.mpl.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.mpl.W, fit2$w)
  expect_identical(fit1$internal$glmdata$.mpl.SUBSET, fit2$subset)

  # Strauss model
  Q <- quadscheme(Xs)
  fit1 <- ppm(Q ~ ., data = exdata, interaction = Strauss(5))
  fit2 <- ppmnet(Q, data = exdata, interaction = Strauss(5))
  expect_identical(as.matrix(fit1$internal$glmdata[,-c(1,2,8)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.mpl.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.mpl.W, fit2$w)
  expect_identical(fit1$internal$glmdata$.mpl.SUBSET, fit2$subset)

  # Geyer model
  Q <- quadscheme(Xg)
  fit1 <- ppm(Q ~ ., data = exdata, interaction = Geyer(5, 1))
  fit2 <- ppmnet(Q, data = exdata, interaction = Geyer(5, 1))
  expect_identical(as.matrix(fit1$internal$glmdata[,-c(1,2,8)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.mpl.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.mpl.W, fit2$w)
  expect_identical(fit1$internal$glmdata$.mpl.SUBSET, fit2$subset)

})

test_that("ppmnet internals are consistent with ppm -- logi", {

  # Poisson model
  Q <- quadscheme.logi(bei)
  fit1 <- ppm(Q ~ ., data = bei.extra, method = "logi")
  fit2 <- ppmnet(Q, data = bei.extra, method = "logi")
  expect_identical(as.matrix(fit1$internal$glmdata[,c(1,2)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.logi.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.logi.w, fit2$w)
  expect_identical(fit1$internal$glmdata$.logi.B, rep.int(fit2$b, n.quad(Q)))
  expect_identical(fit1$internal$glmdata$.logi.ok, fit2$subset)

  # Strauss model
  Q <- quadscheme.logi(Xs)
  fit1 <- ppm(Q ~ ., data = exdata, interaction = Strauss(5), method = "logi")
  fit2 <- ppmnet(Q, data = exdata, interaction = Strauss(5), method = "logi")
  expect_identical(as.matrix(fit1$internal$glmdata[, c(2:5, 1)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.logi.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.logi.w, fit2$w)
  expect_identical(fit1$internal$glmdata$.logi.B, rep.int(fit2$b, n.quad(Q)))
  expect_identical(fit1$internal$glmdata$.logi.ok, fit2$subset)

  # Geyer model
  Q <- quadscheme.logi(Xg)
  fit1 <- ppm(Q ~ ., data = exdata, interaction = Geyer(5, 1), method = "logi")
  fit2 <- ppmnet(Q, data = exdata, interaction = Geyer(5, 1), method = "logi")
  expect_identical(as.matrix(fit1$internal$glmdata[, c(2:5, 1)]), fit2$x)
  expect_identical(fit1$internal$glmdata$.logi.Y, fit2$y)
  expect_identical(fit1$internal$glmdata$.logi.w, fit2$w)
  expect_identical(fit1$internal$glmdata$.logi.B, rep.int(fit2$b, n.quad(Q)))
  expect_identical(fit1$internal$glmdata$.logi.ok, fit2$subset)

})

test_that("ppmnet works - method = 'mpl'", {
  # Poisson model
  Qp <- quadscheme(Xp)
  fit0p <- ppm(Qp ~ ., data = exdata)
  fit1p <- ppmnet(Qp, data = exdata, nlambda = 20, lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0p), coef(fit1p)[, length(fit1p$lambda)],
               tolerance = 0.001)

  # Strauss model
  Qs <- quadscheme(Xs)
  fit0s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5))
  fit1s <- ppmnet(Qs, data = exdata, interaction = Strauss(5),
                  nlambda = 20, lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0s), coef(fit1s)[, length(fit1s$lambda)],
               tolerance = 0.001)

  # Geyer model
  Qg <- quadscheme(Xg)
  fit0g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1))
  fit1g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1),
                  nlambda = 20, lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0g), coef(fit1g)[, length(fit1g$lambda)],
               tolerance = 0.001)

})

test_that("ppmnet works - method = 'logi'", {
  # Poisson model
  Qp <- quadscheme.logi(Xp)
  fit0p <- ppm(Qp ~ ., data = exdata, method = "logi")
  fit1p <- ppmnet(Qp, data = exdata, method = "logi", nlambda = 20,
                  lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0p), coef(fit1p)[, length(fit1p$lambda)],
               tolerance = 0.01)

  # Strauss model
  Qs <- quadscheme.logi(Xs)
  fit0s <- ppm(Qs ~ ., data = exdata, interaction = Strauss(5),
               method = "logi")
  fit1s <- ppmnet(Qs, data = exdata, interaction = Strauss(5),
                  method = "logi", nlambda = 20, lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0s), coef(fit1s)[, length(fit1s$lambda)],
               tolerance = 0.01)

  # Geyer model
  Qg <- quadscheme.logi(Xg)
  fit0g <- ppm(Qg ~ ., data = exdata, interaction = Geyer(5, 1),
               method = "logi")
  fit1g <- ppmnet(Qg, data = exdata, interaction = Geyer(5, 1),
                  method = "logi", nlambda = 20, lambda.min.ratio = 1e-9)
  expect_equal(coef(fit0g), coef(fit1g)[, length(fit1g$lambda)],
               tolerance = 0.01)

})
