
BT <- function (formula = formula(data), data = list(), tweedie.power = 1,
    ABT = TRUE, n.iter = 100, train.fraction = 1, interaction.depth = 4, 
    shrinkage = 1, bag.fraction = 1, colsample.bytree = NULL, 
    keep.data = TRUE, is.verbose = FALSE, cv.folds = 1, folds.id = NULL, 
    n.cores = 1
    , tree.control = rpart.control(xval = 0, maxdepth = (if (!is.null(interaction.depth)) {
        interaction.depth
    } else {
        10
    }), cp = -Inf, minsplit = 2), weights = NULL, seed = NULL, 
    ...) 
{
    if (!is.null(seed)) 
        set.seed(seed)
    the_call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.omit
    mf[[1]] <- as.name("model.frame")
    m <- mf
    mf <- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    respVar <- as.character(attr(Terms, "variables"))[-1][attr(Terms, 
        "response")]
    explVar <- attr(Terms, "term.labels")
    if (!is.null(attr(Terms, "offset"))) {
        stop("Offset are not supported. For Tweedie model with log-link function, weights (=offset)\n         and response rate variable (=Original response variable/Offset) can instead be used.")
    }
    if (is.null(model.weights(mf))) {
        mf$w <- rep(1, nrow(mf))
    }
    else {
        colnames(mf)[names(mf) == "(weights)"] <- "w"
    }
    w <- "w"
    .check_tweedie_power(tweedie.power)
    .check_ABT(ABT)
    .check_n_iter(n.iter)
    .check_train_fraction(train.fraction)
    .check_interaction_depth(interaction.depth)
    .check_shrinkage(shrinkage)
    .check_bag_fraction(bag.fraction)
    .check_colsample_bytree(colsample.bytree, length(explVar))
    .check_keep_data(keep.data)
    .check_is_verbose(is.verbose)
    .check_cv_folds(cv.folds)
    .check_folds_id(folds.id)
    .check_n_cores(n.cores)
    .check_weights(mf$w)
    if (!is.null(interaction.depth) && tree.control$maxdepth != 
        interaction.depth) {
        stop("interaction.depth and maxdepth defined. If interaction.depth is not null it has to be set to maxdepth.")
    }
    setList <- .create_validation_set(mf, train.fraction)
    training.set <- setList$training.set
    validation.set <- setList$validation.set
    rm(setList)
    rm(mf)
    gc()
    if (is.verbose) 
        message("Fit the model on the whole training set. \n")
    BT_full_results <- BT_call(training.set, validation.set, 
        tweedie.power, respVar, w, explVar, ABT, tree.control, 
        train.fraction, interaction.depth, bag.fraction, shrinkage, 
        n.iter, colsample.bytree, keep.data, is.verbose)
    if (!is.null(folds.id)) {
        numFolds <- length(unique(folds.id))
        if (cv.folds != numFolds) 
            warning("CV folds changed from ", cv.folds, " to ", 
                numFolds, " because of levels in folds.id")
        cv.folds <- numFolds
        folds.id <- as.numeric(as.factor(folds.id))
    }
    if (cv.folds == 1) {
        BT_full_results$cv.folds <- cv.folds
        BT_full_results$call <- the_call
        BT_full_results$Terms <- Terms
        BT_full_results$seed <- seed
        return(BT_full_results)
    }
    if (is.verbose) 
        message("Fit the model on the different CV folds. \n")
    folds <- .create_cv_folds(training.set, cv.folds, folds.id, 
        seed)
    if (n.cores > 1) {
        cl <- makeCluster(n.cores)
        clusterExport(cl, varlist = c("training.set", "tweedie.power", 
            "respVar", "w", "explVar", "ABT", "tree.control", 
            "train.fraction", "interaction.depth", "bag.fraction", 
            "shrinkage", "n.iter", "colsample.bytree", "keep.data", 
            "is.verbose"), envir = environment())
        BT_cv_results <- parLapply(cl, seq_len(cv.folds), function(xx) {
            if (!is.null(seed)) 
                set.seed(seed * (xx + 1))
            valIndex <- which(folds == xx)
            trainIndex <- setdiff(1:length(folds), valIndex)
            BT_call(training.set[trainIndex, ], training.set[valIndex, 
                ], tweedie.power, respVar, w, explVar, ABT, tree.control, 
                train.fraction, interaction.depth, bag.fraction, 
                shrinkage, n.iter, colsample.bytree, FALSE, is.verbose)
        })
        on.exit(stopCluster(cl))
    }
    else {
        BT_cv_results <- lapply(seq_len(cv.folds), function(xx) {
            if (!is.null(seed)) 
                set.seed(seed * (xx + 1))
            valIndex <- which(folds == xx)
            trainIndex <- setdiff(1:length(folds), valIndex)
            BT_call(training.set[trainIndex, ], training.set[valIndex, 
                ], tweedie.power, respVar, w, explVar, ABT, tree.control, 
                train.fraction, interaction.depth, bag.fraction, 
                shrinkage, n.iter, colsample.bytree, FALSE, is.verbose)
        })
    }
    class(BT_cv_results) <- "BTCVFit"
    cv_errors <- .BT_cv_errors(BT_cv_results, cv.folds, folds)
    bestIterCV <- which.min(cv_errors)
    predictions <- predict(BT_cv_results, training.set, cv.folds, 
        folds, bestIterCV)
    BT_full_results$cv.folds <- cv.folds
    BT_full_results$folds <- folds
    BT_full_results$call <- the_call
    BT_full_results$Terms <- Terms
    BT_full_results$seed <- seed
    BT_full_results$BTErrors$cv.error <- cv_errors
    BT_full_results$cv.fitted <- predictions
    return(BT_full_results)
}

BT_call <- function (training.set, validation.set, tweedie.power, respVar,
    w, explVar, ABT, tree.control, train.fraction, interaction.depth,
    bag.fraction, shrinkage, n.iter, colsample.bytree, keep.data,
    is.verbose) 
{
    BT <- list()
    if (is.verbose) 
      message("bag.fraction is not used for the initialization fit.")
    init <- BT_callInit(training.set, validation.set, tweedie.power, 
        respVar, w)
    initF <- list(initFit = init$initFit, training.error = init$trainingError, 
        validation.error = init$validationError)
    currTrainScore <- init$currTrainScore
    currValScore <- init$currValScore
    rm(init)
    gc()
    BT <- BT_callBoosting(cbind(training.set, currTrainScore), 
        cbind(validation.set, currValScore), tweedie.power, ABT, 
        tree.control, interaction.depth, bag.fraction, shrinkage, 
        n.iter, colsample.bytree, train.fraction, keep.data, 
        is.verbose, respVar, w, explVar)
    BT$BTInit <- structure(initF, class = "BTInit")
    class(BT) <- "BTFit"
    return(BT)
}

BT_callInit <- function(training.set, validation.set, tweedie.power
    , respVar, w) {
  #overal rate (this suggest response variable is rate, say pi)
  initFit <- sum(training.set[, w] * training.set[, respVar]) /
    sum(training.set[, w])
  # log(common_rate)
  currTrainScore <- rep(log(initFit), nrow(training.set))

  trainingError <- sum(
      BT_devTweedie(training.set[, respVar], exp(currTrainScore)
      , tweedieVal = tweedie.power, w = training.set[, w])
    ) / nrow(training.set)
    currValScore <- NULL
    validationError <- NULL
    if (!is.null(validation.set)) {
        currValScore <- rep(log(initFit), nrow(validation.set))
        validationError <- sum(BT_devTweedie(validation.set[, 
            respVar], exp(currValScore), tweedieVal = tweedie.power, 
            w = validation.set[, w]))/nrow(validation.set)
    }
    return(list(initFit = initFit, currTrainScore = currTrainScore, 
        currValScore = currValScore, trainingError = trainingError, 
        validationError = validationError))
}



BT_callBoosting <- function (training.set, validation.set, tweedie.power, ABT, tree.control, 
    interaction.depth, bag.fraction, shrinkage, n.iter, colsample.bytree, 
    train.fraction, keep.data, is.verbose, respVar, w, explVar){

    sampRow <- 1:nrow(training.set)
    currFormula <- as.formula(paste("residuals ~ ", paste(explVar, collapse = " + ")))
    training.error <- NULL
    validation.error <- NULL
    oob.improvement <- NULL
    listFits <- list()

    for (iTree in seq_len(n.iter)) {
        if (is.verbose) {
            if ((iTree <= 10) || (iTree <= 100 && iTree%%10 == 
                0) || (iTree%%100 == 0)) 
                message("Iteration: ", iTree)
            }
    #reponse varaible (i.e., the obseved rate) devided by exp(current score)
        training.set[, "residuals"] <- training.set[, respVar]/exp(training.set[, "currTrainScore"])
        training.set[, "iWeights"] <- training.set[, w] * (exp(training.set[, "currTrainScore"])^(2 - tweedie.power))
        if (bag.fraction < 1) {
            sampRow <- sample(1:nrow(training.set), bag.fraction * 
                nrow(training.set))
            oobRow <- setdiff(1:nrow(training.set), sampRow)
            oldOOBError <- (sum(BT_devTweedie(training.set[oobRow, 
                respVar], exp(training.set[oobRow, "currTrainScore"]), 
                tweedieVal = tweedie.power, w = training.set[oobRow, 
                  w]))/length(oobRow))
            if (iTree == 1) 
                initOOB <- oldOOBError
        }
        if ((!is.null(colsample.bytree)) && (colsample.bytree != 
            length(explVar))) {
            sampVar <- explVar[sample(1:length(explVar), colsample.bytree)]
            currFormula <- as.formula(paste("residuals ~ ", paste(sampVar, 
                collapse = " + ")))
        }
        if (tweedie.power == 1) {
            currFit <- rpart(currFormula, training.set[sampRow, ]
              , weights = training.set[sampRow, "iWeights"]
              , method = "poisson", control = tree.control, y = FALSE)
        }
        else {
            stop("Currently implemented for Poisson distribution only.")
        }
        if (!is.null(interaction.depth)) {
            if (!ABT) {
                splittingStrategy <- .BT_splittingStrategy(currFit, 
                  interaction.depth)
                if (!is.null(splittingStrategy) && length(splittingStrategy) > 
                  0) 
                  currFit <- snip.rpart(currFit, toss = splittingStrategy)
            }
            else {
                currFit <- prune(currFit, cp = currFit$cptable[, 
                  "CP"][max(which(currFit$cptable[, "nsplit"] <= 
                  interaction.depth))])
            }
        }
        currFit$where <- NULL
        training.set[, "currTrainScore"] <- training.set[, "currTrainScore"] + 
            shrinkage * log(predict(currFit, newdata = training.set, 
                type = "vector"))
        listFits[[iTree]] <- currFit
        training.error[iTree] <- sum(BT_devTweedie(training.set[sampRow, 
            respVar], exp(training.set[sampRow, "currTrainScore"]), 
            tweedieVal = tweedie.power, w = training.set[sampRow, 
                w]))/length(sampRow)
        if (bag.fraction < 1) {
            oob.improvement[iTree] <- (oldOOBError - (sum(BT_devTweedie(training.set[oobRow, 
                respVar], exp(training.set[oobRow, "currTrainScore"]), 
                tweedieVal = tweedie.power, w = training.set[oobRow, 
                  w]))/length(oobRow)))
        }
        if (!is.null(validation.set)) {
            validation.set[, "currValScore"] <- validation.set[, 
                "currValScore"] + shrinkage * log(predict(currFit, 
                newdata = validation.set, type = "vector"))
            validation.error[iTree] <- sum(BT_devTweedie(validation.set[, 
                respVar], exp(validation.set[, "currValScore"]), 
                tweedieVal = tweedie.power, w = validation.set[, 
                  w]))/nrow(validation.set)
        }
    }
    BT_CallBoosting <- list()
    BT_CallBoosting$BTErrors <- structure(list(training.error = training.error, 
        validation.error = validation.error, oob.improvement = oob.improvement), 
        class = "BTErrors")
    class(listFits) <- "BTIndivFits"
    BT_CallBoosting$BTIndivFits <- listFits
    BT_CallBoosting$distribution <- tweedie.power
    BT_CallBoosting$var.names <- explVar
    BT_CallBoosting$response <- respVar
    BT_CallBoosting$w <- w
    if (keep.data) 
        BT_CallBoosting$BTData <- structure(list(training.set = training.set[, 
            !(colnames(training.set) %in% c("iWeights", "residuals"))], 
            validation.set = validation.set), class = "BTData")
    BT_CallBoosting$BTParams <- structure(list(ABT = ABT, train.fraction = train.fraction, 
        shrinkage = shrinkage, interaction.depth = interaction.depth, 
        bag.fraction = bag.fraction, n.iter = n.iter, colsample.bytree = colsample.bytree, 
        tree.control = tree.control), class = "BTParams")
    BT_CallBoosting$keep.data <- keep.data
    BT_CallBoosting$is.verbose <- is.verbose
    BT_CallBoosting$fitted.values <- training.set[, "currTrainScore"]
    return(BT_CallBoosting)
}
<bytecode: 0x10dc1e4d8>
<environment: namespace:BT>