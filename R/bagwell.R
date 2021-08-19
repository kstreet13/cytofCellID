
# this is just a attempted re-creation of the Bagwell method


# source('R/utils.R')
# #fname <- 'data/1_U20S (RKH20201118).fcs'
# fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
# sce <- prepData(fname)


##########################
### find peaks & fit Gaussians in select channels (ordered) 
##########################
# Order (from Bagwell):
# Beads, Offset, Width, Center, DNA1, Residual, Event_length, Live/Dead, DNA2

runBagwell <- function(sce, plot = FALSE){
    # Beads, Offset, Width, Center, DNA1, Residual, Event_length, Live/Dead, DNA2
    
    # set up T stats matrix
    tstats <- matrix(0, nrow = ncol(sce), ncol = 9)
    colnames(tstats) <- c('Beads', 'Offset', 'Width', 'Center', 'DNA1', 'Residual', 'Event_length', 'Live/Dead', 'DNA2')
    keep <- rep(TRUE, ncol(sce))
    
    ### Beads ###
    x <- assay(sce,'exprs')['Bead',]
    if(plot){
        #layout(matrix(1:2, nrow=1))
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Beads - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x)
    fits <- fit_gaussian(x, groups)
    gaussParams <- fits[which.min(fits[,'mean']), ] # closest peak
    tstats[,'Beads'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 1)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Beads - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Offset ###
    x <- log1p(int_colData(sce)$Offset)
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Offset - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'nobs']), ] # largest peak
    tstats[,'Offset'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 2)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Offset - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Width ###
    x <- log1p(int_colData(sce)$Width)
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Width - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'nobs']), ] # largest peak
    tstats[,'Width'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 3)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Width - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Center ###
    x <- log1p(int_colData(sce)$Center)
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Center - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'nobs']), ] # largest peak
    tstats[,'Center'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 4)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Center - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### DNA1 ###
    x <- assay(sce,'exprs')['DNA1',]
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "DNA1 - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    # bagwell: "tallest peak"
    gaussParams <- fits[which.max(fits[,'height']), ] # tallest peak
    tstats[,'DNA1'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 5)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "DNA1 - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Residual ###
    x <- log1p(int_colData(sce)$Residual)
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Residual - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'height']), ] # tallest peak
    tstats[,'Residual'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 6)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Residual - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Event Length ###
    x <- log1p(int_colData(sce)$Event_length)
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Event Length - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'nobs']), ] # largest peak
    tstats[,'Event_length'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 7)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Event Length - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### Viability ###
    x <- assay(sce,'exprs')['Live_Dead',]
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "Viability - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.min(fits[,'mean']), ] # closest peak
    tstats[,'Live/Dead'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 8)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "Viability - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
    }
    
    
    ### DNA2 ###
    x <- assay(sce,'exprs')['DNA2',]
    if(plot){
        d1 <- density(x[which(keep)])
        plot(d1$x, d1$y, type='l', main = "DNA2 - Before", xlab = '', ylab = '')
        polygon(c(min(d1$x),d1$x,max(d1$x)), c(0,d1$y,0), col = 'grey60')
    }
    groups <- find_groups(x[which(keep)])
    fits <- fit_gaussian(x[which(keep)], groups)
    gaussParams <- fits[which.max(fits[,'height']), ] # tallest peak
    tstats[,'DNA2'] <- (x - gaussParams['mean']) / gaussParams['sd']
    keep <- rowSums(.5*tstats^2) < qchisq(.99, 9)
    if(plot){
        d2 <- density(x[which(keep)])
        plot(d2$x, d2$y, xlim = range(d1$x), type='l', main = "DNA2 - After", xlab = '', ylab = '')
        polygon(c(min(d2$x),d2$x,max(d2$x)), c(0,d2$y,0), col = 'grey30')
        layout(1)
    }
    
    return(list(keep = keep, tstats = tstats))
}

# res <- runBagwell(sce)

