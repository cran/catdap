.packageName <- "catdap"

#========================================================================

catdap1 <- function(cdata, response.names = NULL, plot = 1, ask = TRUE) {

#========================================================================

    if (plot != 1 && plot != 2)
      plot <- 0

#   cdata                    # original data
    n <- dim(cdata)[2]       # total number of variable
    nsamp <- dim(cdata)[1]   # sample size
    catdata <- array(0, dim=c(n, nsamp))
    for (i in 1:n)
      catdata[i, ] <- cdata[, i]

#   title                   # variable names
    title <- names(cdata)
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

    dname <- list()
    for (i in 1:n)
      if (is.numeric(cdata[, i]) == FALSE) 
        dname[[i]] <- levels(ordered(cdata[, i]))

    catdap1.out <- catdap01(catdata, title, dname, response.names) 

    if (plot != 0)
      plot.catdap1(catdap1.out, plot, ask)

    class(catdap1.out) <- "catdap1"
    return(catdap1.out)
}

#==========================================================================

catdap1c <- function(ctable, response.names = NULL, plot = 1, ask = TRUE) {

#==========================================================================

    if (plot != 1 && plot != 2)
      plot <- 0

#  ctable                  # original contingency table data
    id <- dim(ctable)
    n <- length(id)
    y <- as.array(ctable)
    yn <- length(y)
    nsamp <- 0
    for (i in 1:yn)
      nsamp <- nsamp + y[[i]]
    cdata <- array(0, dim=c(n, nsamp))

    nn <-  yn
    x <- c(0:(yn-1))
    for (i in n:1) {
      nn <-  nn / id[i]
      nd <- 0
      for (j in 1:yn) {
        m <-  as.integer(x[j] / nn)
        if (y[[j]] != 0) for (k in 1:y[[j]]) cdata[i, nd+k] <- m + 1
        nd <- nd + y[[j]]
        x[j] <- x[j] %% nn
      }
    }

#   title                   # variable names
    title <- names(dimnames(ctable))
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

    dname <- dimnames(ctable)
    if (is.null(dname)) {
      dname <- rep(list(), n)
      for (i in 1:n)
        dname[[i]] <- as.character(1:id[i])
    }

    catdap1c.out <- catdap01(cdata, title, dname, response.names) 

    if (plot != 0)
      plot.catdap1(catdap1c.out, plot, ask)

    class(catdap1c.out) <- "catdap1"
    return(catdap1c.out)
}

#==========================================================

catdap01 <- function(cdata, title, dname, response.names) {

#==========================================================

    n <- dim(cdata)[1]       # total number of variable
    nsamp <- dim(cdata)[2]   # sample size

#   response variable numbers
    ier <- 0
    if (is.null(response.names) ) {
      ires <- c(1:n)
      l <- n
    } else {
      l <- length(response.names)
      ires <- NULL
      m <- 0
      for (i in 1:l)
        for (j in 1:n)
          if (response.names[i] == title[j]) {
            m <- m + 1
            ires <- c(ires, j)
          }
      if (m == 0)
        stop(" Error : response variable name is wrong.")
      if (m < l)
        ier <- 3001
      l <- m
    }

#   conv                    # variable number and recode numbers
    recode <- 0
    iconv <- array(0, dim=c(20, recode))

#   ex                      # variable numbers of explanatory variables
    n1 <- 0
    iex <- 0

#   minmax
    item1 <- rep(0, n)
    item2 <- rep(0, n)
    for (i in 1:n)
      item1[i] <- as.integer(min(cdata[i, ]))
    for (i in 1:n)
      item2[i] <- as.integer(max(cdata[i, ]))
    nc <- rep(0, n)
    for (i in 1:n)
      nc[i] <- item2[i] - item1[i] + 1
    n4 <- max(nc)

    if (length(dname) == 0) {
      for (i in 1:n)
        if (is.numeric(cdata[i, ])) {
          j1 <- item1[i]
          j2 <- item2[i]
          dname[[i]] <- c(j1:j2)
        }
    } else {
      for (i in 1:n)
        if (is.null(dname[[i]]) == TRUE) {
          j1 <- item1[i]
          j2 <- item2[i]
          dname[[i]] <- c(j1:j2)
        }
    }

#   skip                     # variable number and code number to delete
    iskip1 <- 0
    iskip <- 0
    
    z <- .Call("catdap1m",
                 as.integer(nsamp),
                 as.integer(n),
                 as.integer(l),
                 as.integer(recode),
                 as.integer(n1),
                 as.integer(iskip1),
                 as.integer(n4),
                 as.integer(item1),
                 as.integer(item2),
                 as.integer(iconv),
                 as.integer(ires),
                 as.integer(iex),
                 as.integer(iskip),
                 as.integer(cdata))

    nc <- z[[1L]]
    ia <- array(z[[2L]], dim=c(n4, n4, l, n))
    p <- array(z[[3L]], dim=c(n4, n4, l, n))
    total <- array(z[[4L]], dim=c(n4, n))
    aic <- array(z[[5L]], dim=c(l, n))
    for (i in 1:l)  aic[i, ires[i]] <- NA
    ord <- array(z[[6L]], dim=c(l, n-1))
    ia <- aperm(ia, c(2, 1, 3, 4))
    p <- aperm(p, c(2, 1, 3, 4))

    cross.table <- list()
    for (i in 1:l) {
      i1 <- 1
      i2 <- nc[ires[i]]
      cross.table[[i]] <- list()
      for (j in 1:n) {
        if (j == ires[i]){
          cross.table[[i]][[j]] <- NULL
        } else {
          j1 <- 1
          j2 <- nc[j]
          cross.table[[i]][[j]] <- list()
          cross.table[[i]][[j]]$res <- title[ires[i]]
          cross.table[[i]][[j]]$n <- array(ia[j1:j2, i1:i2, i, j], dim=c(j2,i2),
            dimnames = list(dname[[j]], dname[[ires[i]]]))
          cross.table[[i]][[j]]$p <- array(p[j1:j2, i1:i2, i, j], dim=c(j2,i2))
        }
      }  # for (j in 1:n) end
    }  # for (i in 1:l) end

    etotal <- matrix(list(), n)
    for (j in 1:n) {
      j1 <- 1
      j2 <- nc[j]
      etotal[j] <- list(total[j1:j2, j])
    }

    catdap01.out <- list(tway.table = cross.table, title = title, total = etotal,
                         aic = aic, aic.order = ord, ier = ier)

    return(catdap01.out)
}

#=================================

print.catdap1 <- function(x,...) {

#=================================

    l <- dim(x$aic)[1]
    n <- dim(x$aic)[2]
    ires <- rep(0, l)

#--------------------------------------------------------------#
    cat("\n < Summary of AIC's for the two-way tables >\n\n")
#--------------------------------------------------------------#
    cat("         Explanatory \n           variables ")
    for (i in 1:l)
      cat(sprintf("\t%20s", x$tway.table[[i]][[1]]$res))
    cat("\n")
    for (i in 1:(l+2))
      cat("---------------------")
    cat("\n")
    for (jj in 1:n) {
      cat(sprintf("%20s", x$title[jj]))
      for (j in 1:l) {
        if (is.na(x$aic[j,jj]) == TRUE) {
          cat("\t                 -  ")
          ires[j] <- jj
        } else {
          if (is.null(x$aic[j,jj]) == FALSE)
            cat(sprintf("\t            %8.2f", x$aic[j,jj]))
        }
      }
      cat("\n")
    }

    nc <- n-1
    if (nc < 6) {
      nc1 <- nc
      nr <- 1
    } else if (nc < 11) {
      nc1 <- as.integer((nc+1) / 2)
      nr <- 2
    } else {
      nc1 <- 5
      nr <- 2
    }

    dname <- list()
    item <- rep(0, n)
    for (i in 1:n)
      if (i != ires[1]) {
        dname[[i]] <- dimnames(x$tway.table[[1]][[i]]$n)[[1]]
        dname[[ires[1]]] <- dimnames(x$tway.table[[1]][[i]]$n)[[2]]
      }
    for (i in 1:n)
      item[i] <- length(dname[[i]])

    for (i in 1:l) {
      i1 <- 1
      i3 <- ires[i]
      i2 <- item[i3]
#------------------------------------------------------------------------------#
      cat("\n\n < List of explanatory variables arranged in ascending order of AIC >")
#------------------------------------------------------------------------------#
      cat(sprintf("\n\n Response variable  : %s\n\n", x$title[i3]))
      cat("     No.\t  Explanatory\tNumber of categories\t  A I C\t\tDifference of AIC\n")
      cat("        \t  variable   \tof exp. var.\n")
      for (j in 1:9)
        cat("----------")
      cat("\n")
      as <- 0.0
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        aic <- x$aic[i, nv]
        if (j != 1)
          as <- aic - as
        cat(sprintf("%8i %20s\t%8i\t\t%8.2f\t%8.2f\n", j,x$title[nv], item[nv],
            aic, as))
        as <- aic
      }  # for (j in 1:nc) end
#--------------------------------------------------------------------------#
      cat("\n < Two-way tables arranged in ascending order of AIC >\n\n")
#--------------------------------------------------------------------------#

      cat(sprintf("\t\t\t\t( %s )\n", x$title[i3]))
      nsamp <- 0
      ptt <- 0.0
      pt <- rep(0, i2)
      cat("        ")
      for (ii in i1:i2) {
        cat(sprintf("\t%s\t", dname[[i3]][ii]))
        nsamp <- nsamp + x$total[[i3]][ii]
      }
      for (ii in i1:i2) {
        pt[ii] <- x$total[[i3]][ii] / nsamp * 100.0
        ptt <- ptt + pt[ii]
      }

      cat("\t\tTotal\n")
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        j1 <- 1
        j2 <- item[nv]
        cat(sprintf("( %s )\n", x$title[nv]))
        for( k in j1:j2) {
          cat(sprintf("%s\t", dname[[nv]][k]))
          tc <- 0
          for (ii in i1:i2)
            cat(sprintf("%8i ( %5.1f )", x$tway.table[[i]][[nv]]$n[k,ii],
                x$tway.table[[i]][[nv]]$p[k,ii]))
          cat(sprintf("%8i ( %5.1f )\n", x$total[[nv]][k], ptt))
        }
        cat("Total\t")
        for (ii in i1:i2) cat(sprintf("%8i ( %5.1f )", x$total[[i3]][ii], pt[ii]))
        cat(sprintf("%8i ( %5.1f )\n\n", nsamp, ptt))
      }  # for (j in 1:nc) end
    }  # for (i in 1:l) end

    if (x$ier == 3001)  cat(" caution : some response variable name is wrong.\n")

}


#=================================================

plot.catdap1 <- function(x, plot.type, ask, ...) {

#=================================================

    l <- dim(x$aic)[1]
    n <- dim(x$aic)[2]
    ires <- rep(0, l)

    for (jj in 1:n) 
      for (j in 1:l) 
        if (is.na(x$aic[j,jj]) == TRUE)
          ires[j] <- jj

    dname <- list()
    item <- rep(0, n)
    for (i in 1:n)
      if (i != ires[1]) {
        dname[[i]] <- dimnames(x$tway.table[[1]][[i]]$n)[[1]]
        dname[[ires[1]]] <- dimnames(x$tway.table[[1]][[i]]$n)[[2]]
      }

    for (i in 1:n)
      item[i] <- length(dname[[i]])

    old.par <- par(no.readonly=TRUE)
    
    par.ask <- FALSE
    for (i in 1:l) {
      if (plot.type == 1)
        plot.single1(x$title, dname, ires[i], item, x$aic.order[i, ], x$aic[i, ],
                     x$tway.table[[i]], x$total, old.par, par.ask)
      if (plot.type == 2)
        plot.single2(x$title, dname, ires[i], item, x$aic.order[i, ], x$aic[i, ],
                     x$tway.table[[i]], old.par, par.ask)
      par.ask <- ask
    } 

    i1 <- 1
    i3 <- ires[1]
    i2 <- item[i3]
    nsamp <- 0
    for (ii in i1:i2)
      nsamp <- nsamp + x$total[[i3]][ii]
    plot.grshade(x$title, x$aic, ires, nsamp, old.par, ask)
}


#=========================================================================

catdap2 <-
 function(data, pool = NULL, response.name, accuracy = NULL, nvar = NULL,
 additional.output = NULL, missingmark = NULL, pa1 = 1, pa2 = 4, pa3 = 10,
 print.level = 0, plot = 1) {

#=========================================================================

    if (print.level != 1)
      print.level <- 0
    if (plot != 1 && plot != 2)
      plot <- 0

#   data			# original data
    nsamp <- dim(data)[1]	# total number of variable
    n <- dim(data)[2]	# sample size

#   title                   # variable names
    title <- names(data)
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

#   response variable number
    ier <- 0
    ires <- 0
    for (i in 1:n)
      if (response.name == title[i])
        ires <- i
    if (ires == 0)
      stop(" Error : response variable name is wrong.")

#   accuracy                 # accuracy of measurement
    xx <- accuracy
    if (is.null(accuracy))
      xx <- rep(0, n)

#   pooling                    # the way of pooling categories of the variable
    if (is.null(pool))
      pool <- rep(1, n)

#   missing value handling
    typeu <- missingmark
    if (is.null(missingmark))
      typeu <- 0
    z2 <- MissingCount(data, pool, xx, typeu)
    cdata <- z2$cdata
    nmiss <- z2$nmiss
    mix <- z2$mix
    dname <- z2$dname
    pool <- z2$pooling
    xx <- z2$acc
    nsamp <- dim(cdata)[2]
    if (typeu != 0) {
      nmiss.total <- 0
      for (i in 1:length(nmiss))
        nmiss.total <- nmiss.total + sum(nmiss[[i]])
      if (nmiss.total == 0)
        typeu <- 0
    }

#   minmax                   # minimum and maximum of the code values
    item1 <- rep(0, n)
    item2 <- rep(0, n)

    for (i in 1:n) {
      if (i == ires) {
        if (pool[i] < 0) {   # continuous response variable
          item1[i] <- 1
          item2[i] <- -pool[ires]
          if (xx[i] == 0) {
            xx[i] <- (max(cdata[i, ]) - min(cdata[i, ])) / 50
            cat(sprintf("\n Warning : %d-th accuracy is corrected to %8.2f\n",
                i, xx[i]))
          }
          pool[i] <- 0
        } else if (pool[i] < 2) {
          stop(" Error : 'pool' is less than 0 for the response variable.")

        } else {   # categorical response variable
          item1[i] <- as.integer(min(cdata[i, ]))
          item2[i] <- as.integer(max(cdata[i, ]))
          pool[i] <- 2
          xx[i] <- 0
        }
      } else { # if (i != ires)
        if (pool[i] == 2) {
          xx[i] <- 0
          item1[i] <- as.integer(min(cdata[i, ]))
          item2[i] <- as.integer(max(cdata[i, ]))
        } else {
          if (pool[i] < 0)
            pool[i] <- 1
          if (xx[i] == 0) {   # continuous explanatory variable
            xx[i] <- (max(cdata[i, ]) - min(cdata[i, ])) / 50
            cat(sprintf("\n Warning : %d-th accuracy is corrected to %8.2f\n",
                i, xx[i]))
          }
          item1[i] <- 0
          item2[i] <- 0
        } # if(pool[i] == 2) end
      }  # if (i == ires) end
    } # for (i in 1:n) end

    l <- 1
    recode <- 0
    iskip1 <- 0

    it <- 0
    for (i in 1:n)
      if (xx[i] != 0) {
        it <- 1
        break
      }
    xx.in <- rep(0, n)
    for (i in 1:n)
      xx.in[i] <- xx[i]

    iconv <- rep(0, 20)
    face <- ires
    iskip <- rep(0, 20)
    isk <- rep(0, 2)
    sk <- rep(0, 20)
    ida <- array(0, dim=c(n, nsamp))
    da <- array(0, dim=c(n, nsamp))
    if (it == 0)
      ida <- cdata
    if (it == 1)
      da <- cdata

    if (is.null(nvar))
      nvar <- n
    if (nvar <2)
      stop(" Error : 'nvar' is greater than or equal to 2.")
    if (nvar > n)
      stop(" Error : 'nvar' is less than or equal to the number of variable.")
    novv <- 0
    for (i in 1:n)
      if (pool[i] == 0 || pool[i] > 0)
        novv <- novv + 1
    nov <- min(nvar, novv)

#   additional contingency table 
    icl <- 0
    if (is.null(additional.output)) {
      icl <- 0
      icls <- rep(0, max(10, nov+1))
    } else {
      lexp <- length(additional.output)
      icls <- array(0, c(max(10, nov+1), lexp))
      for (i in 1:lexp) {
        nex <- length(additional.output[[i]])
        if (nex > (nov-1)) {
          cat(sprintf("\n Warning : the number of explanatory variable names for additional output is less than nvar %i.\n\n", nov))
        } else if (nex > max((nov-1),8)) {
          cat(sprintf("\n Warning : the number of explanatory variable names for additional output is less than %i.\n\n", max((nov-1), 9)))
        } else {
          icl <- icl + 1
          icls[1,icl] <- nex + 1
          icls[2,icl] <- ires
          m <- 0
          for (j in 1:nex)
            for (k in 1:n)
              if (additional.output[[i]][j] == title[k]) {
                if (k != ires) {
                  icls[2 + j, icl] = k
                  m <- m + 1
                }
              }
            if (m != nex) {
              cat("\n Warning : explanatory variable name for additional output is wrong.\n\n")
              icl <- icl - 1
            }
        }
      }   # for (i in 1:lexp) end
      icls <- icls[, 1:icl]
    }

#   dimension parameters
    n11 <- item2[ires] - item1[ires] + 1
    n2 <- 0
    for (i in 1:n) {
      n1 <- item2[i] - item1[i] + 1
      n2 <- max(n2, n1)
    }
    n33 <- 100
    if (nsamp / 2 + 1 < n33)
      n33 <- nsamp / 2 + 2
    n33 <- max(n33, n2)
    if (it == 0)
      n33=n2

    ikr <- n * pa2
    jkr <- pa3
    icl1 <- icl + 1

### border check
    ikkk <- n33 ** 3 * pa1 * pa1
    n31 <- n * ikkk * icl1
    ddmax <- 2 ** 24 - 1
    ikkkm <- as.integer(ddmax / n / icl1)
    if (ikkk < 2**7)
      ikkk = 2 ** 7
    if (ikkk > ikkkm)
      ikkk <- ikkkm
    
    eps <- 1.0e-10   # an error tolerance in the difference of AIC's

    z <- .Call("catdap2m",
               as.integer(nsamp),
               as.integer(n),
               as.integer(l),
               as.integer(recode),
               as.integer(iskip1),
               as.integer(it),
               as.integer(nov),
               as.integer(icl),
               as.integer(item1),
               as.integer(item2),
               as.integer(pool),
               as.integer(iconv),
               as.integer(ires),
               as.integer(iskip),
               as.integer(isk),
               as.double(sk),
               as.integer(icls),
               as.double(xx.in),
               as.integer(ida),
               as.double(da),
               as.double(typeu),
               as.integer(mix),  
               as.integer(n11),
               as.integer(n33),
               as.integer(ikr),
               as.integer(jkr),
               as.integer(ikkk),
               as.double(eps))

    ier <- z[[22L]]
 
    if (ier[1] == 0) {
      lk77 <- z[[14L]]
      iab <- array(z[[1L]], dim=c(n, n11, n33))

      totalc <- z[[2L]]
      ttrr <- array(z[[3L]], dim=c(n, n33))
      ab <- array(z[[4L]], dim=c(n, n33))
      iaa <- array(z[[5L]], dim=c(n, n33))
      pa <- array(z[[6L]], dim=c(n, n11, n33))
      idata <- z[[7L]]
      ite <- z[[8L]]
      aic <- z[[9L]]
      aaam <- z[[10L]][2:lk77]
      caa <- array(z[[11L]], dim=c(ikr, jkr))
      caa <- caa[1:lk77, 1:jkr]
      icaa <- z[[12L]][1:lk77]
      nnia <- z[[13L]][2:lk77]
      morder <- z[[15L]][1:lk77]

      iby <- array(z[[16L]], dim=c(nov ,ikkk, icl1))
      ibc <- array(z[[17L]], dim=c(n11, ikkk, icl1))
      pbc <- array(z[[18L]], dim=c(n11,ikkk,icl1))
      aic1 <- z[[19L]]
      iabse <- array(z[[20L]], dim=c(nov, n33, icl1))
      baic <- z[[21L]]

      j1 <- 1
      if (is.null(dname)) {
        for (j in 1:n) {
          j2 <- ite[j]
          dname[[j]] <- c(j1:j2)
        }
      } else {
        for (j in 1:n)
          if (is.na(dname[[j]][1])) {
            j2 <- ite[j]
            dname[[j]] <- c(j1:j2)
          }
      }

      iab <- aperm(iab, c(1, 3, 2))
      pa <- aperm(pa, c(1, 3, 2))
      cross.table <- list()
      i1 <- 1
      i2 <- ite[ires]
      for (j in 1:n) {
        if (j == ires){
          cross.table[[j]] <- NA
        } else {
          j1 <- 1
          j2 <- ite[j]
          cross.table[[j]] <- list()
          cross.table[[j]]$res <- title[ires]
          cross.table[[j]]$n <- array(iab[j, j1:j2, i1:i2], dim=c(j2,i2),
                                dimnames = list(dname[[j]], dname[[ires]]))
          cross.table[[j]]$p <- array(pa[j, j1:j2, i1:i2], dim=c(j2,i2))
        }
      }

      total <- matrix(list(), n)
      for (j in 1:n) {
        j1 <- 1
        j2 <- ite[j]
        if (j == ires)
          total[j] <- list(totalc[j1:j2])
        if (j != ires)
          total[j] <- list(ttrr[j, j1:j2])
      }

      nint <- max(ite) + 1
      interval <- list()
      for (j in 1:n) {
        j3 <- idata[j]
        j1 <- 1
        j2 <- ite[j3]
        if (xx[j3] != 0.0)
          interval[[j3]] <- ab[j3, j1:(j2+1)]
        if (xx[j3] == 0.0)
          interval[[j3]] <- as.integer(iaa[j3, j1:j2])
      }
      aic.order <- idata[2:n]

      aorder <- morder- 1
      nsub <- length(aorder)

      nv <- rep(0, nsub)
      ncc <- rep(0, nsub)
      aaic <- rep(0, nsub)
      icaa <- icaa[2:lk77]

      nzero <- 0
      for (i in 1:nsub) {
        j <- aorder[i]
        if (j > 0) {
         nv[i] <- icaa[j] - 1
         ncc[i] <- nnia[j]
         aaic[i] <- aaam[j]
        } else {
          nzero <- nzero + 1
        }
      }
      if (nzero > 1) {
        nsub <- nsub - nzero + 1
        nv <- nv[1:nsub]
        ncc <- ncc[1:nsub]
        aaic <- aaic[1:nsub]
      }

      cexp <- list()
      for (i in 1:nsub) {
        if (aorder[i] > 0) {
          j <- aorder[i] + 1
          k <- nv[i] + 1
          cexp[[i]] <- aic.order[caa[j, 2:k] - 1]
        } else {
          cexp[[i]] <- 0
        }
      }

      iby <- aperm(iby, c(2, 1, 3))
      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      ctbl <- list()
      ctbl.cnum <- list()
      ctbl.p <- list()
      caic <- rep(0, icl1)
      cinterval <- list()
      cexvar <- list()
      idx <- 0

      if (aorder[1] > 0) {
        np1 <- ncc[1]
        nn <- icaa[aorder[1]]
        ctbl.cnum[[icl1]] <- iby[1:np1, 2:nn, 1]
        ctbl[[icl1]] <- ibc[1:np1, i1:i2, 1]
        ctbl.p[[icl1]] <- pbc[1:np1, i1:i2, 1]
      } else {
        ctbl.cnum[[icl1]] <- 0
        ctbl[[icl1]] <- 0
        ctbl.p[[icl1]] <- 0
      }
      caic[icl1] <- aic1[1]

### The output of the additional analysis ### 
      if (icl > 0) {  
        intval <- list()
        for (ic in 1:icl) {
          np1 <- 1
          idx <- idx + 1
          cint <- NULL
          exvar<- NULL
          if (icl == 1)
            nvv <- icls[1] - 1
          if (icl > 1)
            nvv <- icls[1, ic] - 1
          for (i in 1:nvv) {
            if (icl == 1) exv <- icls[i+2]
            if (icl > 1) exv <- icls[i+2, ic]
            intv <- NULL
            j1 <- 1
            j2 <- ite[exv]
            if (xx[exv] != 0.0) {
              p <- iabse[(i+1), j1:(j2+1), 1+ic]
              nint <- length(p)
              for (k in 1:nint) 
                if (p[k] != 0)
                  intv <- c(intv, interval[[exv]][p[k]])
              np1 <- np1 * (length(intv) - 1)
            } else if (xx[exv]==0.0) {
              intv  <- as.integer(iaa[exv, j1:j2])
              np1 <- np1 * j2
            }
            intval[[i]] <- intv
            cint <- c(cint, intval[i])
            exvar <- c(exvar, exv)
          }    # for (i in 1:nvv) end
          cinterval[[idx]] <- cint
          cexvar[[idx]] <- exvar
          if (icl == 1)
            nn <- icls[1]
          if (icl > 1)
            nn <- icls[1, ic]
          ctbl.cnum[[ic]] <- iby[1:np1, 2:nn, 1+ic]
          ctbl[[ic]] <- ibc[1:np1, i1:i2, 1+ic]
          ctbl.p[[ic]] <- pbc[1:np1, i1:i2, 1+ic]
          caic[ic] <- aic1[1+ic]
        }   # for (ic in 1:icl) end
      }  #  if (icl > 0) end


### The output of the Minimum AIC Model ###

      idx <- idx + 1    
      cint <- NULL
      exvar <- NULL

      if (nv[1] > 0) {
        intval <- list()
        for (i in 1:nv[1]) {
          exv <- cexp[[1]][i]
          intv <- NULL
          j1 <- 1
          j2 <- ite[exv]
          if (xx[exv] != 0.0) {
            p <- iabse[(i+1), j1:(j2+1), 1]
            nint <- length(p)
            for (k in 1:nint) 
              if (p[k] != 0)
                intv <- c(intv, interval[[exv]][p[k]])
            np1 <- np1 * (length(intv) - 1)
          } else if (xx[exv]==0.0) {
            intv  <- as.integer(iaa[exv, j1:j2])
            np1 <- np1 * j2
          }
          intval[[i]] <- intv
          cint <- c(cint, intval[i])
          exvar <- c(exvar, exv)
        }    # for (i in 1:nv[1]) end
      }
      cinterval[[idx]] <- cint
      cexvar[[idx]] <- exvar

      iaddflg <- 0
      nbest <- 1
      if (nsub > 1) { 
        aicm1 <- aic1[1]
        for (i in 2:nsub) {
          if (iaddflg == 0) { 
           aicm2 <- aaic[i]
            daic <- abs(aicm2 - aicm1) / max(abs(aicm2), abs(aicm1))
            if (daic > eps)
              break
            nbest <- nbest + 1

#          if (iaddflg == 0) {
            for (j in 1:n)
              xx.in[j] <- xx[j]

            z1 <- addaicm(i, nsamp, n, l, recode, iskip1, it, nov, item1, item2,
                  pool, iconv, ires, iskip, isk, sk, xx.in, ida, da, typeu, mix,
                  n11, n33, pa1, ikr, jkr, ikkk, ikkkm, cexp, eps)

            if (z1$ier == 0) {
              iaddflg <- 1
              idx <- idx + 1
              np1 <- ncc[nbest]
              idflg <- is.null(dim(z1$ctable$cnum[[1]]))
              if (idflg == TRUE)
                ccnum <- z1$ctable$cnum[[1]][1:np1]
              if (idflg == FALSE)
                ccnum <- z1$ctable$cnum[[1]][1:np1, ]
              ctbl.cnum[[idx]] <- ccnum
              ctbl[[idx]] <- z1$ctable$n[[1]]
              ctbl.p[[idx]] <- z1$ctable$p[[1]]
              caic <- c(caic, aicm2)
              iabse <- z1$iabse

              cint <- NULL
              intval <- list()
              exvar <- NULL
              for (j in 1:nv[i]) {
                exv <- cexp[[i]][j]
                intv <- NULL
                j1 <- 1
                j2 <- ite[exv]
                if (xx[exv] != 0.0) {
                  p <- iabse[(j+1), j1:(j2+1)]
                  nint <- length(p)
                  for (k in 1:nint) 
                    if (p[k] != 0 )
                      intv <- c(intv, interval[[exv]][p[k]])
                } else if (xx[exv]==0.0) {
                  intv  <- as.integer(iaa[exv, j1:j2])
                }
                intval[[j]] <- intv
                cint <- c(cint, intval[j])
                exvar <- c(exvar, exv)
              } # for (j in 1:nv[i]) end
              cinterval[[idx]] <- cint
              cexvar[[idx]] <- exvar
            } # if (z1$ier == 0) end
          } # if (iaddflg == 0) end 
        }  # if (i in 2:nsub) end
      }  # if (nsub > 2) end

      catdap2.out <- list(title = title, accuracy = xx, ires = ires,
                      print.level = print.level, tway.table = cross.table,
                      total = total, interval = interval, base.aic = baic,
                      aic = aic, aic.order = aic.order, nsub = nsub, nv = nv,
                      ncc = ncc, aaic = aaic, cexp = cexp,
                      ctable = list(n = ctbl, p = ctbl.p, cnum = ctbl.cnum),
                      ctable.interval = list(exvar = cexvar, range = cinterval),
                      caic = caic, missing = nmiss)

      if (plot != 0)
        plot.catdap2(catdap2.out, plot)    

    } else { # if (ier[1] != 0)
      if (ier[1] == 2048) { # pa1
        kkj <- ier[2]
        if (kkj < ikkkm) {
          pa1n <- ceiling( sqrt(kkj / (n33 ** 3)) )
          ier[2] <- pa1n
          ikkkn <- n33 ** 3 * pa1n * pa1n
          if (ikkkn > ikkkm) ier[2] <- -999    
        } else {
          ier[2] <- -999
        }  # if (kkj < ikkkm) end
      }  # if (ier[1] == 2048) end
      catdap2.out <- list(ier=ier)
    }

    class(catdap2.out) <- "catdap2"
    return(catdap2.out)

}


#================================================================

convi <- function(data) {

#================================================================

    x1 <- unique(data)
    y1 <- rank(x1)
    n1 <- length(data)
    n2 <- length(y1)
    cdata <- rep(0, n1)
    for (i in 1:n1)
      for (j in 1:n2)
        if (data[i] == x1[j])
          cdata[i] <- y1[j]
    cname <- sort(x1)
    return(list(cdata = cdata, cname = cname))
} 

#================================================================

addaicm <- 
 function(iex2, nsamp, n, l, recode, iskip1, it, nov, item1, item2, pool,
 iconv, ires, iskip, isk, sk, xx, ida, da, typeu, nmiss, n11, n33, pa1, ikr,
 jkr, ikkk, ikkkm, cexp, eps) {

#================================================================

    icl <- 1
    icl1 <- icl + 1
    nex <- length(cexp[[iex2]])
    icls <- c(nex+1, ires)
    for (i in 1:nex)
      icls <- c(icls, cexp[[iex2]][i])

    z <- .Call("catdap2m",
               as.integer(nsamp),
               as.integer(n),
               as.integer(l),
               as.integer(recode),
               as.integer(iskip1),
               as.integer(it),
               as.integer(nov),
               as.integer(icl),
               as.integer(item1),
               as.integer(item2),
               as.integer(pool),
               as.integer(iconv),
               as.integer(ires),
               as.integer(iskip),
               as.integer(isk),
               as.double(sk),
               as.integer(icls),
               as.double(xx),
               as.integer(ida),
               as.double(da),
               as.double(typeu),
               as.integer(nmiss),  
               as.integer(n11),
               as.integer(n33),
               as.integer(ikr),
               as.integer(jkr),
               as.integer(ikkk),
               as.double(eps))

    ier <- z[[22L]]
 
    if (ier[1] == 0 ) {
      ite <- z[[8L]]
      iby <- array(z[[16L]], dim=c(nov ,ikkk, icl1))
      ibc <- array(z[[17L]], dim=c(n11, ikkk, icl1))
      pbc <- array(z[[18L]], dim=c(n11,ikkk,icl1))
      iabse <- array(z[[20L]], dim=c(nov, n33, icl1))

      iby <- aperm(iby, c(2, 1, 3))
      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      i1 <- 1
      i2 <- ite[ires]

      np1 <- 1
      for (i in 1:(icls[1]-1))
        np1 <- np1 * ite[icls[2+i]]
      nn <- icls[1]

      jby <- list(iby[1:np1, 2:nn, 2])
      jbc <- list(ibc[1:np1, i1:i2, 2])
      qbc <- list(pbc[1:np1, i1:i2, 2])
      iabse <- iabse[, , 2]

      return(list(ctable = list(n = jbc, p = qbc, cnum = jby),
                  iabse = iabse, ier = 0))

    } else {
      if (ier[1] == 2048) { # pa1
        kkj <- ier[2]
        if (kkj < ikkkm) {
          pa1n <- ceiling( sqrt(kkj / (n33 ** 3)))
          ikkkn <- n33 ** 3 * pa1n * pa1n
          ier[2] <- pa1n
          if (ikkkn > ikkkm)
            ier[2] <- -999    
        } else {
          ier[2] <- -999
        }  # if (kkj < ikkkm) end
        if (ier[2] == -999)
          cat("Error : Working area for second contingency table is too short. pa1 can no longer set a larger value.\n")
        if ( ier[2] != -999)
          cat(sprintf(" Error : Working area for second contingency table is too short, try pa1= %i.\n", pa1n))
      }  # if (ier[1] == 2048) end
      return(list(ier = ier[1]))
    }
}

#==================================

print.catdap2 <- function(x, ...) {

#==================================

    eps <- 1.0e-10   # an error tolerance in the difference of AIC's
    dl10 <- "----------"
    bl7 <- "       "
    bl14 <- "              "

    if (is.null(x$ier[1])) {
      n <- length(x$title)
      res <- x$ires
      lk <- x$nsub
      ktt <- max(x$nv)

      dname <- list()
      nc <- rep(0, n)
      for (i in 1:n)
        if (i != res) { 
          dname[[i]] <- dimnames(x$tway.table[[i]]$n)[[1]]
          dname[[res]] <- dimnames(x$tway.table[[i]]$n)[[2]]
        }
      for (i in 1:n)
        nc[i] <- length(dname[[i]])

#----------------------------------------------------------------------------#
      cat("\n<< List of single explanatory variables (arranged in ascending order of AIC) >>\n")
#----------------------------------------------------------------------------#
      cat(sprintf(" Response variable : %s \t(base AIC = %8.2f)\n", x$title[res],
          x$base.aic))
      for (nrep in 1:8)
        cat(dl10)
      cat("\n")
      cat(bl14, bl14, "Number of\n")
      cat(bl7, "      Explanatory     categories", bl14, " Difference\n")
      cat(bl7, "      variables       of exp. var.    A I C      of AIC",
          bl7, "Weight\n")
      for (nrep in 1:8)
        cat(dl10)
      cat("\n")

      nrank <- 0
      daic2 <- 1.0e+5

      for (i in 1:(n-1)) {
        idx <- x$aic.order[i]
        aaaa <- x$aic[idx]
        if (i == 1)
          aaa1 <- x$aic[idx]
        daic <- aaaa - aaa1
        w <- exp(-1. / 2 * daic)

        if (i > 1) {
          if (aaaa==0.0e0 && aaa2==0.0e0) {
            daic2 <- 0.0e0
          } else {
            daic2 <- abs(aaaa - aaa2) / max(abs(aaaa), abs(aaa2))
          }
        }

        if (daic2 > eps)
          nrank <- i
        cat(sprintf("%5i", nrank))
        cat(sprintf("%20s", x$title[idx]))
        cat(sprintf("     %8i     %8.2f     %8.2f     %8.2f\n", nc[idx],
            x$aic[idx], daic, w))
        aaa2 <- aaaa
      }

#-----------------------------------------------------------------------#
      cat("\n\n<< Two-way tables arranged in ascending order of AIC >>\n\n")
#-----------------------------------------------------------------------#

      i <- 1
      i1 <- 1
      i2 <- nc[res]
      nsamp <- sum(x$total[[res]][i1:i2])
      ptc <- rep(0, i2)
      for (i in i1:i2)
        ptc[i] <- (x$total[[res]][i] * 100.) / nsamp
      ptt <- sum(ptc[i1:i2])

      nmiss <- x$missing[[res]]
      ntype <- length(nmiss)
      if (ntype==1 && nmiss[1]==0)
        ntype <- 0
      if (x$accuracy[res] != 0) 
        print.Note(x$title[res], x$accuracy[res], x$interval[[res]], dname[[res]],
                   i1, i2, ntype, nmiss)

      cat(sprintf("                     ( %15s )\n", x$title[res]))
      for (j in 1:(n-1)){
        j3 <- x$aic.order[j]
        j1 <- 1
        j2 <- nc[j3]
        nmiss <- x$missing[[j3]]
        ntype <- length(nmiss)
        if (ntype==1 && nmiss[1]==0)
          ntype <- 0

        cat(sprintf("                 ", x$title[j3]))
        for (ii in i1:i2)
          cat(sprintf("%15s   ", dname[[res]][ii]))

        cat("          Total\n")
        cat(sprintf("( %15s )\n", x$title[j3]))
        for (ii in i1:(i2+2))
          cat("-----------------")
        cat("-----\n")
        for (k in j1:j2) {
          ptr <- 0
          if (x$total[[j3]][k] != 0)
            ptr <- (x$total[[j3]][k] * 100.) / x$total[[j3]][k]
          cat(sprintf("%8i           ", k))
          for (ii in i1:i2)
            cat(sprintf("%8i ( %5.1f )", x$tway.table[[j3]]$n[k,ii],
                x$tway.table[[j3]]$p[k,ii]))
          cat(sprintf("%8i ( %5.1f )\n", x$total[[j3]][k], ptr))
        }  # for (k in j1:j2) end

        for (ii in i1:(i2+2))
        cat("-----------------")
        cat("-----\n     Total         ")
        for (ii in i1:i2)
          cat(sprintf("%8i ( %5.1f )", x$total[[res]][ii], ptc[ii]))
        cat(sprintf("%8i ( %5.1f )\n\n", nsamp, ptt))

        nmiss <- x$missing[[j3]]
        ntype <- length(nmiss)
        if (ntype==1 && nmiss[1]==0)
          ntype <- 0
        print.Note(x$title[j3], x$accuracy[j3], x$interval[[j3]], dname[[j3]],
                   j1, j2, ntype, nmiss)

      }  # for (j in 1:(n-1)) end

      if (x$print.level == 0) {
#-----------------------------------------------------------------------------#
        cat("\n<< AIC's of the models with k explanatory variables (k=1,2,...) >>\n")
#-----------------------------------------------------------------------------#
        for (j in 1:ktt) {
          ijk <- 0
          nrank <- 0
          daic2 <- 1.0e+5

          for (ij in 1:lk) {
            lk4 <- x$nv[ij]
            if (lk4 == j) {
              aaaa <- x$aaic[ij]
              if (ijk == 0)
                aaic1 <- x$aaic[ij]
              daic <- aaaa - aaic1
              w <- exp(-1. / 2 * daic)
              ijk <- ijk + 1

              if (ijk == 1) {
                cat(sprintf("\n Number of explanatory variables = %i\n", j))
                for (nrep in 1:8)
                  cat(dl10)
                cat("\n")
                cat(bl14, bl14, "Number of\n")
                cat(bl7, "      Explanatory     categories", bl14,
                    " Difference\n")
                cat(bl7,
                    "      variables       of exp. var.    A I C      of AIC",
                    bl7, "Weight\n")
                for (nrep in 1:8)
                  cat(dl10)
                cat("\n")
              }

              if (ijk > 1) {
                if (aaaa==0.0e0 && aaic2==0.0e0) {
                  daic2 <- 0.0e0
                } else {
                  daic2 <- abs(aaaa - aaic2) / max(abs(aaaa), abs(aaic2))
                }
              }

              if (daic2 > eps)
                nrank <- ijk
              if (lk4 == 1 || lk4 > 1)  {
                cat(sprintf("%5i", nrank))
                lkk <- x$cexp[[ij]][1]
                cat(sprintf("%20s", x$title[lkk]))
                cat(sprintf("     %8i", x$ncc[ij]))
                cat(sprintf("     %8.2f", x$aaic[ij]))
                cat(sprintf("     %8.2f", daic))
                cat(sprintf("     %8.2f\n", w))
              }

              if (lk4 == 2 || lk4 > 2)
                for (i in 2:lk4) {
                  lkk <- x$cexp[[ij]][i]
                  cat(sprintf("     %20s\n", x$title[lkk]))
                }

            }  # if (lk4 == j) end
            aaic2 <- aaaa
          }  # for (ij in 1:lk) end
        }  # for (j in 1:ktt) end
      }  # if (x$print.level == 0) end

#-------------------------------------------------------------------#
      cat("\n\n<< Summary of subsets of explanatory variables >>\n\n")
#-------------------------------------------------------------------#

      cat(sprintf(" Response variable : %s\n\n", x$title[res]))
      for (nrep in 1:8)
        cat(dl10)
      cat("\n")
      cat(bl14, bl14, "Number of\n")
      cat(bl7, "      Explanatory     categories", bl14, " Difference\n")
      cat(bl7,
          "      variables       of exp. var.    A I C      of AIC",
          bl7, "Weight\n")
      for (nrep in 1:8)
        cat(dl10)
      cat("\n")
      ijk <- 0
      nrank <- 0
      daic2 <- 1.0e+5
      aaic2 <- -1.0e+5
      if (x$print.level == 0)
       lsub <- lk
      if (x$print.level == 1)
        lsub <- min(lk, 30)
      for (ij in 1:lsub) {
        lk4 <- x$nv[ij]
        if (lk4 > 0) {
          aaaa <- x$aaic[ij]
          if (ijk == 0 )
            aaic1 <- x$aaic[ij]
          daic <- aaaa - aaic1
          w <- exp(-1. / 2 * daic)
          ijk <- ijk + 1
          if (ijk > 1) {
            if (aaaa==0.0e0 && aaic2==0.0e0) {
              daic2 <- 0.0e0
            } else {
              daic2 <- abs(aaaa - aaic2) / max(abs(aaaa), abs(aaic2))
            }
          }

          if (daic2 > eps)
            nrank <- ijk
          cat(sprintf("%5i", nrank))
          lkk <- x$cexp[[ij]][1]
          cat(sprintf("%20s", x$title[lkk]))
          cat(sprintf("     %8i", x$ncc[ij]))
          cat(sprintf("     %8.2f", x$aaic[ij]))
          cat(sprintf("     %8.2f", daic))
          cat(sprintf("     %8.2f\n", w))

          if (lk4 == 2 || lk4 > 2)
            for (i in 2:lk4) {
              lkk <- x$cexp[[ij]][i]
              cat(sprintf("     %20s\n", x$title[lkk]))
            }
          aaic2 <- aaaa

        } else if (aaic2 < 0) {
          ijk <- ijk + 1
          cat(sprintf("%5i", ijk))
          cat("               - - -")
          aaaa <- 0
          if (ijk == 1 ) {
            aaic1 <- 0
            daic <- 0
            w <- 1.0
          } else {
            daic <- aaaa - aaic1
            w <- exp(-1. / 2 * daic)
          }
          cat(sprintf("     %8i", 0))
          cat(sprintf("     %8.2f", 0.0))
          cat(sprintf("     %8.2f", daic))
          cat(sprintf("     %8.2f\n", w))

          aaic2 <- aaaa
        }  # if (lk4 > 0) end
      }  # for (ij in 1:lsub) end

      nctbl <- length(x$caic)
      icl <- nctbl - 1
      if ((icl != 0) && (x$caic[nctbl] == x$caic[nctbl-1]))
        icl <- icl - 1
      for (k in 1:nctbl) {
        lk4 <- 0
        if (length(x$ctable.interval$exvar) > (k-1))
          lk4 <- length(x$ctable.interval$exvar[[k]])
        if (k > icl) {

#------------------------------------------------------------------------------#
          cat("\n\n<< Contingency table constructed by the best subset of explanatory variables >>\n\n")
#------------------------------------------------------------------------------#
          if (lk4 == 0)
            icflg <- 0
          idm <- 0
        } else {
#------------------------------------------------------------------------------#
          cat("\n\n<< The output of the additional analysis >>\n\n")
#------------------------------------------------------------------------------#
        }

        cat(sprintf(" X(1) : %s\n", x$title[res]))
        etitle <- list()
        if (lk4 > 0) {
          j1 <- 1
          j2 <- rep(0, lk4)
          mc <- rep(0, lk4)
          for (i in 1:lk4) {
            lkk <- x$ctable.interval$exvar[[k]][i]
            cat(sprintf(" X(%i) : %s\n", (i+1), x$title[lkk]))
            etitle <- c(etitle, x$title[lkk])
            j2[i] <- length(x$ctable.interval$range[[k]][[i]])
            if (x$accuracy[lkk] > 0)
              j2[i] <- j2[i] - 1
          }
          cat("\n")

          mp <- 1
          for (i in 1:lk4)
            mp <- mp * j2[i]
          idm <- j2[1]
          if (lk4 > 1)
            for (i in 2:lk4)
              idm <- c(idm, j2[i])
          for (i in 1:lk4)
            cat(" X  ")
        }   # if (lk4 > 0) end

        cat("\t\t response variable X(1)\n")

        if (lk4 > 0)
          for (i in 1:lk4)
            cat(sprintf("(%i) ", (i+1)))

        for (i in i1:i2)
          cat(sprintf("           %3i      ", i))
        cat("     Total\n")
        for (ii in i1:(i2+2))
          cat("-----------------")
        cat("\n")

        idf <- TRUE
        if (lk4> 0) {
          np <- dim(x$ctable$cnum[[k]])[1]
          if (is.null(np)) np <- length(x$ctable$cnum[[k]])

          idf <- is.null(dim(x$ctable$cnum[[k]])) 
          for (i in 1:np) {
            if (idf == TRUE)
              cat(sprintf(" %2i ", x$ctable$cnum[[k]][i]))
            if (idf == FALSE)
              for (j in 1:lk4)
                cat(sprintf(" %2i ", x$ctable$cnum[[k]][i,j]))
            for (j in i1:i2)
              cat(sprintf("%8i ( %5.1f )", x$ctable$n[[k]][i,j],
                  x$ctable$p[[k]][i,j]))
            ibct <- sum(x$ctable$n[[k]][i, i1:i2])
            pbct <- sum(x$ctable$p[[k]][i, i1:i2])
            cat(sprintf("%8i ( %5.1f )\n", ibct, pbct))
          }  # end for (i in 1:np)
          for (ii in i1:(i2+2))  cat("-----------------")
        }   # if (lk4 > 0) end

        cat("\n   Total")
        if (idf == FALSE && lk4 > 2)
          for (j in 1:(lk4-2))
            cat("    ")
        for (j in i1:i2 )
          cat(sprintf("%8i ( %5.1f )", x$total[[res]][j], ptc[j]))
        cat(sprintf("%8i ( %5.1f )\n", nsamp, ptt))

        nmiss <- x$missing[[res]]
        ntype <- length(nmiss)
        if (ntype==1 && nmiss[1]==0)
          ntype <- 0
        print.CtableNote(1, x$title[[res]], x$accuracy[[res]], x$interval[[res]],
                         dname[[res]], i1, i2, ntype, nmiss)

        if (lk4 > 0) {
          for (i in 1:lk4) {
            lkk <- x$ctable.interval$exvar[[k]][i]

            nmiss <- x$missing[[lkk]]
            ntype <- length(nmiss)
            if (ntype==1 && nmiss[1]==0)
              ntype <- 0

            ix <- i + 1
            print.CtableNote(ix, x$title[[lkk]], x$accuracy[[lkk]],
                             x$ctable.interval$range[[k]][[i]], dname[[lkk]], j1,
                             j2[i], ntype, nmiss)

          }   #  for (i in 1:lk4) end
        }   # if (lk4 > 0) end

        cat(sprintf("\nAIC = %8.2f\n", x$caic[k]))
        cat(sprintf("base AIC = %8.2f\n\n", x$base.aic))

      } # for (k in  1:ntbl) end

      nbest <- 1
      if (lk > 1) { 
        aicm1 <- x$aaic[1]
        for (i in 2:lk) {
          aicm2 <- x$aaic[i]
          daic <- abs(aicm2 - aicm1) / max(abs(aicm2), abs(aicm1))
          if (daic > eps)
            break
          nbest <- nbest + 1
        }
        if (nbest == 3)
          cat("<NOTE> There is another subset with the minimum AIC.\n\n")
        if (nbest > 3)
          cat(sprintf("<NOTE> There are another %3i subsets with the minimum AIC.\n\n", nbest-2))
      }
    } else {
        print.err2(x$ier)
    }

}


#============================================

plot.catdap2 <- function(x, plot.type, ...) {

#============================================

    n <- length(x$title)
    res <- x$ires
    lk <- x$nsub

    dname <- list()
    nc <- rep(0, n)
    for (i in 1:n)
      if (i != res) { 
        dname[[i]] <- dimnames(x$tway.table[[i]]$n)[[1]]
        dname[[res]] <- dimnames(x$tway.table[[i]]$n)[[2]]
      }
    for (i in 1:n)
        nc[i] <- length(dname[[i]])

    i1 <- 1
    i2 <- nc[res]

    nctbl <- length(x$caic)
    icl <- nctbl - 1
    if ((icl != 0) && (x$caic[nctbl] == x$caic[nctbl-1]))
      icl <- icl - 1

    icflg <- 1
    for (k in 1:nctbl) {
      ctable.dname <- list()
      if (x$accuracy[[res]] == 0.0)
        ctable.dname[[1]] <- dname[[res]]
      if (x$accuracy[[res]] != 0.0)
        ctable.dname[[1]] <- c(i1:i2)

      lk4 <- 0
      if (length(x$ctable.interval$exvar) > (k-1))
        lk4 <- length(x$ctable.interval$exvar[[k]])
      if (k > icl) {
        if (lk4 == 0)
          icflg <- 0
        idm <- 0
      }

      etitle <- list()
      if (lk4 > 0) {
        j1 <- 1
        j2 <- rep(0, lk4)
        for (i in 1:lk4) {
          ix <- i + 1
          lkk <- x$ctable.interval$exvar[[k]][i]
          etitle <- c(etitle, x$title[lkk])
          j2[i] <- length(x$ctable.interval$range[[k]][[i]])
          if (x$accuracy[lkk] > 0)
            j2[i] <- j2[i] - 1
          if (x$accuracy[[lkk]] == 0.0)
            ctable.dname[[ix]] <- dname[[lkk]]
          if (x$accuracy[[lkk]] != 0.0)
            ctable.dname[[ix]] <- c(j1:j2[i])
        } 
        mp <- 1
        for (i in 1:lk4)
          mp <- mp * j2[i]
        idm <- j2[1]
        if (lk4 > 1)
          for (i in 2:lk4)
            idm <- c(idm, j2[i])
      }

      plot.mosaic(plot.type, res, icl, k, mp, idm, x$title, etitle, dname,
                  x$ctable$n, ctable.dname, nc, x$aic.order, x$aic,
                  x$tway.table, x$total, icflg)

    } # for (k in  1:ntbl) end


}


#==============================================================================

print.Note <- function(title, accuracy, interval, dname, i1, i2, ntype, nmiss) {

#==============================================================================

      cat("<Note>\n")
      cat(sprintf(" %s\n", title))
      if (accuracy != 0.0) {
        cat("\tcategory    \tvalue range\n")
        if (ntype == 0) {
          for (i in i1:i2)
            cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", i, interval[i],
                interval[i+1]))
        } else if (ntype != 0) {
          for (i in i1:(i2-ntype))
            cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", i, interval[i],
                interval[i+1]))
          for (i in 1:ntype)
            cat(sprintf("\t%8i    \tmissing of type %i\n", i2-ntype+i, nmiss[i]))
        }
    } else {
        cat("\tcategory    \tvariable value\n")
        if (ntype == 0) {
          for (i in i1:i2)
            cat(sprintf("\t%8i    \t%s\n", i, dname[i]))
        } else if (ntype != 0) {
          for (i in i1:(i2-ntype))
            cat(sprintf("\t%8i    \t%s\n", i, dname[i]))
          for (i in 1:ntype)
            cat(sprintf("\t%8i    \tmissing of type %i\n", i2-ntype+i, nmiss[i]))
      }
    }
    cat("\n\n")
}

#=================================================================

print.CtableNote <- function(ix, title, accuracy, interval, dname, i1, i2,
                             ntype, nmiss) {

#=================================================================

    if (ix == 1 )
      cat("\n<Note>\n")
    cat(sprintf("X(%i) : %s\n", ix, title))

    if (accuracy != 0.0) {
      cat("\tcategory    \tvalue range\n")
      if (ntype == 0) {
        for (i in i1:i2) 
          cat(sprintf("\t%8i    \t%12.5e   -   %12.5e\n", i, interval[i],
              interval[i+1]))
      } else if (ntype != 0) {
        for (i in i1:(i2-ntype))
          cat(sprintf("\t%8i    \t%12.5e   -   %12.5e\n", i, interval[i],
              interval[i+1]))
        for (i in 1:ntype)
          cat(sprintf("\t%8i    \tmissing of type %i\n", i2-ntype+i, nmiss[i]))
      }
    } else {
      cat("\tcategory    \tvariable value\n")
      if (ntype == 0) {
        for (i in i1:i2)
          cat(sprintf("\t%8i    \t%s\n", i, dname[i]))
      } else if (ntype != 0) {
        for (i in i1:(i2-ntype))
          cat(sprintf("\t%8i    \t%s\n", i, dname[i]))
        for (i in 1:ntype)
          cat(sprintf("\t%8i    \tmissing of type %i\n", i2-ntype+i, nmiss[i]))
      }
    }
}


#============================

print.err2 <- function(err) {

#============================

    ier <- err[1]
    eval <- err[2]

    if (ier == 2002) {
      cat(" Error : Working area for multi-dimensional table is not enough. Try larger pa2.\n")
    } else if (ier == 2003) {
      cat(" Error : Working area for multi-dimensional table is not enough. Try larger pa3.\n")
    } else if (ier == 2048) {
      if (eval != -999)
        cat(sprintf(" Error : Working area for contingency table is too short, try pa1= %i.\n", eval))
      if (eval == -999)
        cat("Error : Working area for contingency table is too short. pa1 can no longer set a larger value.\n")
    } else if (ier == 2037) {
      cat(" caution : the program catdap cannot deal with data sets where the number of non-zero frequency categories of the response variables is less than 2.\n")
    } else if (ier == 2035) {
      cat(" Error : the value of variable is beyond the interval specified in 'min' and 'max'.\n")
    } else if (ier == 2588) {
      cat(" Error :  lk5 > n-1 \n")
    } else if (ier == 650) {
      cat(" Error : the value of 'nvar' is too small for the additional analysis.\n")
    }
}


#======================================================================

# plot function (plot.single1, plot.single2, plot.mosaic, plot.grshade)

#======================================================================

plot.single1 <- function(title, dname, res, item, aic.order, aic, tway.table,
                         total, old.par, ask) {

    n <- length(title)
    i2 <- item[res]

    if (n < 6) {
      nc <- n
      nr <- 1
    } else if (n < 11) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- 5
      nr <- 2
    }
    m <- nc * nr

    new.mai <- old.par$mai
    new.mai[1] <- new.mai[1] * 0.65
    new.mai[3] <- new.mai[3] * 0.5
    newcex.main <- old.par$cex.main * 0.8  
    mtitle <- "Single Explanatory Models\nin ascending order of AIC"

    nplot <- 1
    irflag <- 0
    nw <- 0
    if (ask ==TRUE)
      par(ask=TRUE)

    j1 <- 1
    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tway.table[[iex]]$n

      if (aic[iex] > 0 && irflag == 0) {
        if (nplot%%m == 1) {
          if (j > 1 || nw > 0)
            if (.Device != "null device" || ask == FALSE)
              dev.new()
          par(mfcol=c(nc, nr), mai = new.mai)
          nw <- nw + 1
        }
        irflag <- 1
        y <- total[[res]]
        y <- as.array(y)
        xlabel <- paste(title[res], " Total ", total[[res]][1])
        if (i2 > 1)
          for (ii in 2:i2)
            xlabel <- paste(xlabel, ":", total[[res]][ii])
        dimnames(y) <- list(dname[[res]])
        mosaicplot(y, color = TRUE, main = "", ylab = "AIC= 0.0\n ",
                   xlab = xlabel, off = 0)
        nplot <- nplot + 1
      }

      if (nplot%%m == 1) {
        if (j > 1 || nw > 0)
          if (.Device != "null device" || ask == FALSE)
            dev.new()
        par(mfcol=c(nc, nr), mai = new.mai, cex.main = newcex.main)
        nw <- nw + 1
      }

      j2 <- item[iex]
      y <- array(0, dim=c(j2, i2))
      ncc <- 0
      for (k in j1:j2) {
          ncc <- ncc + 1
          for (ii in 1:i2)
            y[ncc, ii] <- ctbl[k, ii]
        }
      dimnames(y)[1] <- list(dimnames(ctbl)[[1]])
      dimnames(y)[2] <- list(dname[[res]])
      aaic <- round(aic[iex], 2)
      mosaicplot(y, color = TRUE, main = mtitle,
                 ylab = paste("AIC=", aaic, "\n", title[iex]),
                 xlab = paste(title[res]),  dir = c("h", "v"), off = c(10, 0))
      mtitle <- ""
      nplot <- nplot + 1
    }  # for (j in 1:(n-1)) end

    if (irflag == 0) {
      y <- total[[res]]
      y <- as.array(y)
      xlabel <- paste(title[res], " Total ", total[[res]][1])
      if (i2 > 1)
        for (ii in 2:i2)
          xlabel <- paste(xlabel, ":", total[[res]][ii])
      dimnames(y) <- list(dname[[res]])
      mosaicplot(y, color = TRUE, main = "", ylab = "AIC= 0.0\n ",
                 xlab = xlabel, off = 0)
    }

    par(old.par)

}


plot.single2 <- function(title, dname, res, item, aic.order, aic, tway.table,
                         old.par, ask) {

    n <- length(title)
    i2 <- item[res]

    n1 <- n - 1
    if (n1 < 6) {
      nc <- n1
      nr <- 1
    } else if (n1 < 11) {
      nc <- as.integer(n / 2)
      nr <- 2
    } else {
      nc <- 5
      nr <- 2
    }
    m <- nc * nr

    new.mai <- old.par$mai
    new.mai[1] <- new.mai[1] * 0.65
    new.mai[3] <- new.mai[3] * 0.5
    newcex.main <- old.par$cex.main * 0.8  
    mtitle <- "Single Explanatory Models\nin ascending order of AIC"

    nplot <- 1
    nw <- 0
    if (ask ==TRUE)
      par(ask=TRUE)

    j1 <- 1
    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tway.table[[iex]]$n

      if (nplot %% m == 1) {
        if (j > 1 || nw > 0)
        if (.Device != "null device" || ask == FALSE)
          dev.new()
        par(mfcol=c(nc, nr), mai = new.mai, cex.main = newcex.main)
        nw <- nw + 1
      }

      j2 <- item[iex]
      y <- array(0, dim=c(i2, j2))
      for (ii in 1:i2) 
        for (k in j1:j2) y[ii, k] <- ctbl[k, ii]
      dimnames(y)[1] <- list(dimnames(ctbl)[[2]])
      dimnames(y)[2] <- list(dimnames(ctbl)[[1]])
      aaic <- round(aic[iex], 2)
      mosaicplot(y, color = TRUE, main = mtitle,
                 ylab = paste("AIC=", aaic, "\n", title[iex]),
                 xlab = paste(title[res]))
      mtitle <- ""
      nplot <- nplot + 1
    } 

    par(old.par)

}


plot.mosaic <- function(iplt, res, icl, k, mp, idm, title, etitle, dname,
                        ctable, ctable.dname, nc, aic.order, aic, tway.table,
                        total, iflg)
{
  old.par <- par(no.readonly = TRUE)
  ask <- FALSE

  if (k == 1) {
    if (iplt == 1)
      plot.single1(title, dname, res, nc, aic.order, aic, tway.table, total,
                   old.par, ask)
    if (iplt == 2)
      plot.single2(title, dname, res, nc, aic.order, aic, tway.table, old.par,
                   ask)
  }

  if (iflg > 0 ) {
    i1 <- 1
    i2 <- nc[res]
    xtitle <- title[res]
    ytitle <- etitle[1]

    nex <- length(idm)
    jdm <- array(0, nex)
    for (i in 1:nex)
      jdm[i] <- idm[nex-i+1]
    idm <- c(jdm, i2)
    tt <- ctable[[k]]
    t <- array(tt, dim=idm)

    if (iplt == 1) {
      t <- aperm(t, c(nex:1, nex+1))
      for (i in 1:nex)
        dimnames(t)[i] <- list(ctable.dname[[i+1]])
      dimnames(t)[nex+1] <- list(ctable.dname[[1]])

      split <- "h"
      offset <- 10
      if (nex > 1)
        for (i in 2:nex) {
          ytitle <- c(ytitle, etitle[i])
          split <- c(split, "h")
          offset <- c(offset, 10)
        }
      split <- c(split, "v")
      offset <- c(offset, 0)

    } else if (iplt == 2) {
      t <- aperm(t, c(nex+1, nex:1))
      dimnames(t) <- ctable.dname

      if (nex > 1)
        for (i in 2:nex) {
          if (i %% 2 == 0)
            xtitle <- c(xtitle, etitle[i])
          if (i %% 2 == 1)
            ytitle <- c(ytitle, etitle[i])
        }

      split <- NULL
      offset <- NULL

    } # else if (iplt == 2) end

    if (icl > 0 && k < icl+1)  {
      mtitle <- "Mosaicplot for additinal contingency table"
    } else {
      mtitle <- "Minimum AIC Model of\n the Response Variable Distribution"
    }

    dev.new()
    par(mfcol=c(1, 1))
    mosaicplot(t, color = TRUE, main = mtitle, xlab = paste(xtitle),
               ylab = paste(ytitle), dir = split, off = offset)
  } # if (ifg > 0) end

    par(old.par)
}

plot.grshade <- function(title, aic, ires, nsamp, old.par, ask) {

  new.mai <- old.par$mai
  new.mai[1] <- new.mai[1] * 1.5
  new.mai[2] <- new.mai[2] * 1.5
  new.mai[3] <- new.mai[3] * 2.0
  new.mai[4] <- new.mai[4] * 1.0

  response <- NULL
  explanatory <- NULL
  nres <- length(ires)
  for (i in nres:1)
    response <- c(response, title[ires[i]])
  nex <- length(title)
  for (i in 1:nex)
    explanatory <- c(explanatory, title[i])

  d1 <- dim(aic)[1]
  d2 <- dim(aic)[2]
  aics <- aic / nsamp
  aicsmin <- 0.0
  aicsmax <- 0.0
  for (i1 in 1:d1)
    for (i2 in 1:d2) 
      if (is.na(aics[i1, i2]) == FALSE) {
        aicsmin <- min(aicsmin, aics[i1, i2])
        aicsmax <- max(aicsmax, aics[i1, i2])
      }
  if (aicsmin  > -0.1)
    aicsmin <- -0.15
  if (aicsmax < 0)
    aicsmax <- 0.5

  aics1 <- array(0, c(d1, d2))
  for (i in d1:1)
    aics1[i, ] <- aics[d1-i+1, ]
  aics1 <- aperm(aics1, c(2, 1))
  aics2 <- array(0, c(d2+1, d1))
  aics2[1:d2, ] <- aics1[1:d2, ]
  x <- (1:nrow(aics2))
  y <- (1:ncol(aics2))

  bp <- c(aicsmin, -0.1, -0.05, -0.01, 0, aicsmax) - 1.0e-10

  mtitle1 <- "Gray Shading Display of All the AIC's\n\nGray level : "
  mtitle2 <- "N=Number of observations\n 5 \ : \ \ \ \ \ \ \ \ \ \ \ \ \ "
  mtitle3 <- "AIC/N < -0.10\ \n 4 \ : -0.10 =< AIC/N < -0.05\ \n 3 \ : -0.05 =<"
  mtitle4 <- " AIC/N < -0.01\ \n 2 \ : -0.01 =< AIC/N < \ 0.0\ \ \n 1 \ :\ \ "
  mtitle5 <- "0.0\ \ \ =< AIC/N\ \ \ \ \ \ \ \ \ \ \ \ \n"
  main.title <- paste(mtitle1, mtitle2, mtitle3, mtitle4, mtitle5, sep="")

  if (ask == TRUE) {
    par(ask = TRUE, mai = new.mai, cex.axis = 0.8)
  } else {
    dev.new()
    par(mai = new.mai, cex.axis = 0.8)
  }

  image(x, y, aics2, col=gray(0:4/4), breaks=bp, axes = FALSE, xlab="", ylab="",
        main=main.title, sub="explanatory variable", cex.main=0.8, cex.sub=0.8)

  legend("topright", c("5","4","3","2","1"), pch = 15, col = gray(0:4/4),
         bg = "snow2", cex = 0.8, pt.cex = 1.1)

  axis(1, 1:d2, explanatory, las = 3)
  axis(2, 1:d1, response, las = 1)
  box()

  par(old.par)
}


#=====================================================================

Barplot2WayTable <- function (vname, resvar, exvar=NULL, tway.table,
                              interval=NULL) {

#=====================================================================

  nv <- length(vname)
  if (is.null(exvar) == TRUE)
    exvar <- vname
  nresv <- length(resvar)
  nexv <- length(exvar)
  nc <- min(nexv, nv-1)

  if (length(tway.table) != nv) { #  catdap1
    ic <- 1
  } else if (length(tway.table[[1]]) == nv) {  # catdap1
    ic <- 1
  } else {  # catdap2
    ic <- 2
    if (is.null(interval) == TRUE)
      stop(" Error : Class interval is not specified.")
  }

  if (ic == 1) {  # catdap1() or catdap1c()
    nw <- 0
    for (j in 1:nresv) {
      res <- 0
      for (i in 1:nv) 
        if(vname[i] == resvar[j]) {
          res <- i
          break
        }
      if (res == 0)
        stop(" Error : response variable name is wrong.")

      nbar <- 0
      if (nexv > 1)
        par(mfcol = c(nc, 1))
      for (k in 1:nexv) {
        ex <- 0
        for (i in 1:nv)
          if(vname[i] == exvar[k]) {
            ex <- i
            break
          }
        if (ex == 0)
          stop(" Error : explanatory variable name is wrong.")

        if (res != ex) {
          if (nbar == 0 && nw > 0) {
            dev.new()
            if (nexv > 1)
              par(mfcol = c(nc, 1))
          }
          h <- tway.table[[j]][[ex]]$n
          barplot(h, names.arg = dimnames(h)[[2]], legend = dimnames(h)[[1]],
                  args.legend = list(x="topright", title=vname[ex], cex=0.8),
                  space = 0, xlab = vname[res])
          nbar <- nbar + 1
        }
      }  # for (k in 1:nexv) end
      nw <- nw + 1
    }  # for (j in 1:nresv) end

  } else {  # catdap2()
    res <- 0
    for (i in 1:nv)
      if(vname[i] == resvar) {
        res <- i
        break
      }
    if (res == 0)
      stop(" Error : response variable name is wrong.")

    if (nexv > 1)
      par(mfcol = c(nc, 1))
    for (j in 1:nexv) {
      ex <- 0
      for (i in 1:nv)
        if(vname[i] == exvar[j]) {
          ex <- i
          break
        }
      if (ex == 0)
        stop(" Error : explanatory variable name is wrong.")

      if (res != ex) {
        h <- tway.table[[ex]]$n
        nr <- length(dimnames(h)[[1]])
        if (j == 1) {
          nrres <- length(dimnames(h)[[2]])
          resint <- dimnames(h)[[2]]
          resname <- tway.table[[ex]]$res
        }
        nint <- length(interval[[ex]])
        if (nr == nint) {
          barplot(h, names.arg = resint, legend = dimnames(h)[[1]],
                  args.legend = list(x="topright", title=vname[ex], cex=0.8),
                  space = 0, xlab = vname[res])
        } else {
          itext <- NULL
          for (k in 1:nr)
            itext[[k]] <- paste(interval[[ex]][k], "-", interval[[ex]][k+1])
          barplot(h, names.arg = resint, legend = itext,
                  args.legend = list(x="topright", title=vname[ex], cex=0.8),
                  space = 0, xlab = vname[res])
        }
      }
    }  # for (j in 1:nexv) end

    nint <- length(interval[[res]])
    if (nint > nrres) {
      cat("\n<Note>\n")
      cat(sprintf(" %s\n", resname))
      cat("\tcategory    \tvalue range\n")
      for (k in 1:nrres)
        cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", k, interval[[res]][k],
            interval[[res]][k+1]))
    }
  }
  par(mfcol = c(1,1))
}


#=====================================================================

MissingCount <- function (data, pool, xx, typeu) {

#=====================================================================

  nsamp <- dim(data)[1]	# total number of variable
  n <- dim(data)[2]	# sample size

  cdata <- array(0, dim=c(n, nsamp))
  dname <- list()
  nmiss <- list()
  mix <- rep(0, n)
  n1 <- 0

  for (i in 1:n) {
    z <- MissingCount1(data[,i], pool[i], xx[i], typeu, n1)
    pool[i] <- z$pooling
    xx[i] <- z$acc
    cdata[i,] <- z$cdata
    dname[[i]] <- z$dname
    nmiss[[i]] <- z$nmiss
    mix[i] <- length(nmiss[[i]])
    if (mix[i]==1 && nmiss[[i]][1]==0)
      mix[i] <- 0 
    n1 <- z$n1
  }

  if (n1 == n)
    dname <- NULL

    return(list(pooling = pool, acc = xx, cdata = cdata, nmiss = nmiss,
                mix = mix, dname = dname))

}

#=====================================================================

MissingCount1 <- function (data, pool, xx, typeu, n1) {

#=====================================================================

  nsamp <- length(data)
  cdata <- rep(0, nsamp)
  dname <- NA
  nmiss <- 0

#
#     numerical characters
#
  if (is.numeric(data)) {
    if (pool != 2) {
      n1 <- n1 + 1
      cdata <- data
      if (typeu != 0) {
        for (j in 1:nsamp)
          if (cdata[j] >= typeu) {
            k <- as.integer(cdata[j] / typeu)
            nmiss <- c(nmiss, k)
          }
        nmiss <- unique(nmiss)
      }
      if (length(nmiss) > 1)
        nmiss <- nmiss[-1]
      sord <- order(nmiss)
      nmiss <- nmiss[sord]
#
#     categorical numerical characters
#
    } else if (pool == 2) {
      tmp <- data
      if (typeu != 0) {
        for (j in 1:nsamp)
          if (tmp[j] >= typeu) {
            k <- as.integer(tmp[j] / typeu)
            tmp[j] <- as.integer(k * typeu)
            nmiss <- c(nmiss, k)
          }
        nmiss <- unique(nmiss)
      }
      cv <- convi(tmp)
      cdata <- cv$cdata
      dname <- cv$cname
      if (length(nmiss) > 1) {
        nmiss <- nmiss[-1]
        sord <- order(nmiss)
        nmiss <- nmiss[sord]
        nv <- length(dname)
        nm <- length(nmiss)
        for (j in 1:nm)
          dname[nv-j+1] <- paste("missing of type",
                                 as.integer(dname[nv-j+1])/typeu)
      }
    }
#
#     categorical letters
#
  } else {
    tmp <- data
    cdata <- as.numeric(tmp)
    pool <- 2
    xx <- 0.0
    dname <- levels(tmp)
  }

  return(list(pooling = pool, acc = xx, cdata = cdata, nmiss = nmiss,
              dname = dname, n1 = n1))

}

