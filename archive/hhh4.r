library(surveillance)

#forest_fit <- readRDS("data/forest_fit_prop_gb200.rds")
#colnames(forest_fit)
#dt <- filter(forest_fit, replication == 1)

#dt <- arrange(dt, i, year)
#n  <- with(dt, length(unique(i)))

#obs <- matrix(dt$y0, ncol = n)
#colnames(obs) <- unique(dt$i)
#rownames(obs) <- min(dt$year):max(dt$year)


dt_sts <- sts(
  observed = obs
  , start = c(min(dt$year), 1)
  , frequency = 1
  , population = matrix(dt$y0_hat, ncol = n)
)
dim(dt_sts)

model_1 <- list(
  end = list(f = ~ ri() - 1
    , offset = population(dt_sts) #baseline rate
  )
  , ar = list(f = ~  - 1) #occasional outbreaks
  , ne = list(f = ~ - 1)
  , family = "Poisson"
)


#start <- Sys.time()
#set.seed(0203)
#fit1 <- hhh4(dt_sts, model_1)
#end <- Sys.time()
#print(end - start)
