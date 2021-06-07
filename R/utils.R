# used by Bagwell, us
find_groups <- function(x){
    d <- density(x, adjust = 2)
    # initial peaks:
    #  (1) higher than it's 10 neighbors
    #  (2) at leat 1% of the maximum peak height
    peakTF <- sapply(1:length(d$y), function(i){
        neighbors <- c((i-5):(i-1), (i+1):(i+5))
        neighbors <- neighbors[neighbors %in% 1:length(d$y)]
        
        return(all(d$y[i] > d$y[neighbors]) &
                   d$y[i] > .01*max(d$y))
    })
    peaks <- which(peakTF)
    # find valleys (lowpoints between peaks)
    if(sum(peakTF) > 1){
        valleys <- sapply(1:(sum(peakTF)-1), function(i){
            p1 <- which(peakTF)[i]
            p2 <- which(peakTF)[i+1]
            btwn <- p1:p2
            return(btwn[which.min(d$y[btwn])])
        })
    }
    # peak filtering
    # remove peak if not >120% of both nearby valleys
    for(ii in 1:length(peaks)){
        ht <- d$y[peaks[ii]]
        if(ii == 1){
            lv <- 1
        }else{
            lv <- valleys[ii-1]
        }
        if(ii == length(peaks)){
            rv <- length(d$x)
        }else{
            rv <- valleys[ii]
        }
        lv <- d$y[lv]
        rv <- d$y[rv]
        if(ht < 1.2*lv | ht < 1.2*rv){
            peakTF[peaks[ii]] <- FALSE
            print(paste('removing',ii))
        }
    }
    peaks <- which(peakTF)
    # re-find valleys
    if(sum(peakTF) > 1){
        valleys <- sapply(1:(sum(peakTF)-1), function(i){
            p1 <- which(peakTF)[i]
            p2 <- which(peakTF)[i+1]
            btwn <- p1:p2
            return(btwn[which.min(d$y[btwn])])
        })
        out <- cut(x, breaks = c(-Inf, d$x[valleys], Inf))
    }else{
        out <- factor(rep(1, length(x)))
    }
    levels(out) <- 1:length(peaks)
    return(out)
}

# used by Bagwell
fit_gaussian <- function(x, groups){
    d <- density(x)
    t(sapply(levels(groups), function(lv){
        rng <- range(x[which(groups == lv)])
        group.m <- mean(x[which(groups == lv)])
        group.sd <- sd(x[which(groups == lv)])
        group.n <- sum(groups == lv)
        group.ht <- max(d$y[which(d$x >= rng[1] & d$x <= rng[2])])
        return(c(mean = group.m, 
                 sd = group.sd*sqrt(group.n / (group.n-1)), 
                 nobs = group.n, height = group.ht))
    }))
}
