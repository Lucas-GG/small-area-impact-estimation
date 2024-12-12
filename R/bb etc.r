bb <- \(n) {
  U <- runif(n - 1)
  G <- diff(c(0, sort(U), 1))
  G
}

bb(10)
att_ln <- \(m, dt = dat) {
  y1 <- predict(m, mutate(dt, trt = 1), type = "response")
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  mean(log(y1 / y0)[dt$trt == 1])
}

att <- \(m, dt = dat) {
  # "response" is the defualt for fixest, but just to be explicit
  y1 <- predict(m, mutate(dt, trt = 1), type = "response")
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  mean((y1 - y0)[dt$trt == 1])
}

att_risk <- \(m, dt = dat) {
  # "response" is the defualt for fixest, but just to be explicit
  y1 <- predict(m, mutate(dt, trt = 1), type = "response")
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  sum((y1 - y0)[dt$trt == 1]) / sum(dt$n[dt$trt == 1]) * 10^5
}

att_arisk <- \(m, dt = dat) {
  # "response" is the defualt for fixest, but just to be explicit
  y1 <- predict(m, mutate(dt, trt = 1), type = "response")
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  mean(((y1 - y0) / dt$n)[dt$trt == 1]) * 10^5
}

att0 <- \(m, dt = dat) {
  # "response" is the defualt for fixest, but just to be explicit
  y1 <- dt$Y
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  mean((y1 - y0)[dt$trt == 1])
}

att0_risk <- \(m, dt = dat) {
  # "response" is the defualt for fixest, but just to be explicit
  y1 <- dt$Y
  y0 <- predict(m, mutate(dt, trt = 0), type = "response")
  mean(((y1 - y0) / dt$n)[dt$trt == 1]) * 10^5
}
