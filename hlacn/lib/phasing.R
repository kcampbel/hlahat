mufreq.bayes.mt <- function(mufreq, depth.ratio, cellularity,
    ploidy, avg.depth.ratio, weight.mufreq = 100, weight.ratio = 100,
    CNt.min = 1, CNt.max = 7, CNn = 2,
    priors.table = data.frame(CN = CNt.min:CNt.max, value = 1),
    priors.mt = 1) {
    mufreq.tab <- data.frame(F = mufreq, ratio = depth.ratio,
        weight.mufreq = weight.mufreq, weight.ratio = weight.ratio)
    mufreq_types <- mufreq.types.matrix(CNt.min = CNt.min,
        CNt.max = CNt.max, CNn = CNn)
    model.pts <- mufreq.model.points(cellularity = cellularity,
        ploidy = ploidy, mufreq_types = mufreq_types,
        avg.depth.ratio = avg.depth.ratio)
    model.pts <- cbind(mufreq_types, model.pts)
    rows.x <- 1:nrow(mufreq.tab)
    priors <- rep(1, nrow(model.pts))
    for (i in 1:nrow(priors.table)) {
        priors[model.pts$CNt ==
            priors.table$CN[i]] <- priors.table$value[i]
    }
    priors <- priors * priors.mt
    priors <- priors / sum(priors)
    bayes.fit <- function (x, mat, model.pts, priors) {
        test.ratio <- model.pts$depth.ratio
        test.mufrq <- model.pts$mufreqs
        min.offset <- 1e-323
        score.r <- depth.ratio.dbinom(size = mat[x, ]$weight.ratio,
            depth.ratio = mat[x, ]$ratio, test.ratio)
        score.m <- mufreq.dbinom(mufreq = mat[x, ]$F,
            depth.t = mat[x, ]$weight.mufreq, test.mufrq)
        score.r <- score.r * priors
        score.m <- score.m
        post.model <- score.r * score.m
        post.model[post.model == 0] <- min.offset
        res.cn <- model.pts$CNt[which.max(score.r)]
        idx.pts <- model.pts$CNt == res.cn
        model.lik <- cbind(model.pts[idx.pts, 1:3], log(post.model[idx.pts]))
        if (is.null(dim(model.lik))) {
            max.post <- model.lik
        } else {
            max.post <- model.lik[which.max(model.lik[,4]),]
        }
        max.post
    }
    print(cbind(model.pts, priors))
    types.L <- mapply(FUN = bayes.fit, rows.x, MoreArgs = list(
        mat = mufreq.tab, model.pts = model.pts,
        priors = priors), SIMPLIFY = FALSE)
   types.L <- do.call(rbind, types.L)
   types.L

   colnames(types.L) <- c("CNn", "CNt", "Mt", "LPP")
   types.L
}

prior_creator_mt <- function(prior.mt, CNt.min, CNt.max, CNn = 2) {
    # Mt Prior Creator
    # Generates prior weights for combinations of Mt and CNt
    # Prior.mt is a data frame where each row is a combination of Mt, CNt, and desired prior
    # weight. For example, a row with Mt=1, CNt=2, prior.weight=10 will generate a vector of 
    # weights where solutions with these values are weighted 10 times higher than other solutions.
    # Args:
    #     prior.mt(df): Mt, CNt, prior.weight
    #     CNt.min(int): Tumor total copy number min
    #     CNt.max(int): Tumor total copy number max
    #     CNn(int): Normal total copy number
    #
    # Returns:
    #     Vector of prior weights for each possible Mt x CNt solution 
    mufreq_types <- mufreq.types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNn = CNn)
    prior.out <- rep(1, nrow(mufreq_types))
    for(ii in 1:nrow(prior.mt)){
        Mt <- prior.mt[ii,]$Mt
        prior.weight <- prior.mt[ii,]$prior.weight
        CNt <- prior.mt[ii,]$CNt
        prior.out[
            which(mufreq_types$Mt == Mt & mufreq_types$CNt == CNt)
            ] <- prior.weight
    }
    return(prior.out)
}

mufreq_params <- function(
        mufreq,
        depth.ratio,
        cellularity,
        ploidy = ploidy,
        avg.depth.ratio = 1,
        weight.mufreq = 100,
        weight.ratio = 100,
        CNt.min = 5,
        CNt.max = 5,
        CNn = 2,
        priors.table = data.frame(CN = CNt.min:CNt.max, value = 1),
        priors.mt = c(1,2,1,1,2,1)
    )
    {
    p = list(
        mufreq = mufreq,
        depth.ratio = depth.ratio,
        cellularity = cellularity,
        ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio,
        weight.mufreq = weight.mufreq,
        weight.ratio = weight.ratio,
        CNt.min = CNt.min,
        CNt.max = CNt.max,
        CNn = CNn,
        priors.table = priors.table,
        priors.mt = priors.mt
    )
    return(p)
}
