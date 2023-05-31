#=========================================================================

catdap2 <-
 function(data, pool = NULL, response.name, accuracy = NULL, nvar = NULL,
 additional.output = NULL, missingmark = NULL, pa1 = 1, pa2 = 4, pa3 = 10,
 print.level = 0, plot = 1, gray.shade = FALSE) {

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
      stop(" Response variable name is wrong.")

#   accuracy                 # accuracy of measurement
    xx <- accuracy
    if (is.null(accuracy))
      xx <- rep(0, n)

#   pooling                    # the way of pooling categories of the variable
    if (is.null(pool))
      pool <- rep(1, n)

#   missing value handling
    if (is.null(missingmark)) {
      typeu <- 0
    } else {
      if (missingmark < 1) {
        missingmark <- NULL
        typeu <- 0
      } else {
        typeu <- missingmark
      }
    }
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
            warning(gettextf(" %d-th accuracy is corrected to %8.2f", i, xx[i]),
                    domain = NA)
          }
          pool[i] <- 0
        } else if (pool[i] == 2) {   # categorical response variable
          item1[i] <- as.integer(min(cdata[i, ]))
          item2[i] <- as.integer(max(cdata[i, ]))
          pool[i] <- 2
          xx[i] <- 0
        } else {
          stop(" 'pool' is 2 or less than 0 for the response variable.")
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
            warning(gettextf(" %d-th accuracy is corrected to %8.2f", i, xx[i]),
                    domain = NA)
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
      stop(" 'nvar' is greater than or equal to 2.")
    if (nvar > n)
      stop(" 'nvar' is less than or equal to the number of variable.")
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
          warning(gettextf("the number of explanatory variable names for additional output is less than nvar %i", nov), domain = NA)
        } else if (nex > max((nov-1),8)) {
          warning(gettextf("the number of explanatory variable names for additional output is less than %i", max((nov-1), 9)), domain = NA)
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
              warning("explanatory variable name for additional output is wrong")
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

    z <- .Fortran(C_catdap2m,
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
                  iab = integer(n*n11*n33),
                  totalc = integer(n11),
                  ttrr = integer(n*n33),
                  ab = double(n*n33),
                  iaa = integer(n*n33),
                  pa = double(n*n11*n33),
                  idata = integer(n),
                  ite = integer(n),
                  dx = double(n),
                  aaam = double(ikr),
                  caa = integer(ikr*jkr),
                  icaa = integer(ikr),
                  nnia = integer(ikr),
                  lk77 = integer(1),
                  morder = integer(ikr),
                  ibc = integer(n11*ikkk*icl1),
                  pbc = double(n11*ikkk*icl1),
                  aic1 = double(icl1),
                  iabse = integer(nov*n33*icl1),
                  baic = double(1),
                  as.integer(n11),
                  as.integer(n33),
                  as.integer(ikr),
                  as.integer(jkr),
                  as.integer(ikkk),
                  as.double(eps),
                  nrange = integer(nov*icl1),
                  ier = integer(3))

    ier <- z$ier
 
    if (ier[1] == 0) {
      lk77 <- z$lk77
      iab <- array(z$iab, dim=c(n, n11, n33))
      totalc <- z$totalc
      ttrr <- array(z$ttrr, dim=c(n, n33))
      ab <- array(z$ab, dim=c(n, n33))
      iaa <- array(z$iaa, dim=c(n, n33))
      pa <- array(z$pa, dim=c(n, n11, n33))
      idata <- z$idata
      ite <- z$ite
      aic <- z$dx
      aaam <- z$aaam[2:lk77]
      caa <- array(z$caa, dim=c(ikr, jkr))
      caa <- caa[1:lk77, 1:jkr]
      icaa <- z$icaa[1:lk77]
      nnia <- z$nnia[2:lk77]
      morder <- z$morder[1:lk77]

      ibc <- array(z$ibc, dim=c(n11, ikkk, icl1))
      pbc <- array(z$pbc, dim=c(n11, ikkk, icl1))
      aic1 <- z$aic1
      iabse <- array(z$iabse, dim=c(nov, n33, icl1))
      baic <- z$baic
      nr <- array(z$nrange, dim=c(nov, icl1))

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
          cross.table[[j]]$exv <- title[j]
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
      vnames <- list()
      for (i in 1:nsub) {
        if (aorder[i] > 0) {
          j <- aorder[i] + 1
          k <- nv[i] + 1
          cexp[[i]] <- aic.order[caa[j, 2:k] - 1]
          vnames[[i]] <- title[cexp[[i]]]
        } else {
          cexp[[i]] <- 0
          vnames[[i]] <- NULL
        }
      }

      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      ctbl <- list()
      ctbl.p <- list()
      caic <- rep(0, icl1)
      cinterval <- list()
      cexvar <- list()
      nrange <- NULL
      idx <- 0

      if (aorder[1] > 0) {
        np1 <- ncc[1]
        nn <- icaa[aorder[1]]
        ctbl[[icl1]] <- ibc[1:np1, i1:i2, 1]
        ctbl.p[[icl1]] <- pbc[1:np1, i1:i2, 1]
      } else {
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
          exvar <- NULL
          exrange <- NULL
          if (icl == 1)
            nvv <- icls[1] - 1
          if (icl > 1)
            nvv <- icls[1, ic] - 1
          for (i in 1:nvv) {
            if (icl == 1) exv <- icls[i+2]
            if (icl > 1) exv <- icls[i+2, ic]
            exnr <- nr[i, ic+1]
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
            exrange <- c(exrange, exnr)
          }    # for (i in 1:nvv) end
          cinterval[[idx]] <- cint
          cexvar[[idx]] <- exvar
          nrange[[idx]] <- exrange
          if (icl == 1)
            nn <- icls[1]
          if (icl > 1)
            nn <- icls[1, ic]
          ctbl[[ic]] <- ibc[1:np1, i1:i2, 1+ic]
          ctbl.p[[ic]] <- pbc[1:np1, i1:i2, 1+ic]
          caic[ic] <- aic1[1+ic]
        }   # for (ic in 1:icl) end
      }  #  if (icl > 0) end


### The output of the Minimum AIC Model ###

      idx <- idx + 1    
      cint <- NULL
      exvar <- NULL
      exrange <- NULL

      if (nv[1] > 0) {
        intval <- list()
        for (i in 1:nv[1]) {
          exv <- cexp[[1]][i]
          exnr <- nr[i, 1]
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
          exrange <- c(exrange, exnr)
        }    # for (i in 1:nv[1]) end
      }
      cinterval[[idx]] <- cint
      cexvar[[idx]] <- exvar
      nrange[[idx]] <- exrange

      iaddflg <- 0
      nbest <- 1
      if (nsub > 1) { 
        aicm1 <- aic1[1]
        for (i in 2:nsub) {
          if (iaddflg == 0) {
            aicm2 <- aaic[i]
#    2020/09/06
            if (aicm1 == 0 && aicm2 == 0)
              break

            daic <- abs(aicm2 - aicm1) / max(abs(aicm2), abs(aicm1))
            if (daic > eps)
              break
            nbest <- nbest + 1

            for (j in 1:n)
              xx.in[j] <- xx[j]

            z1 <- addaicm(i, nsamp, n, l, recode, iskip1, it, nov, item1, item2,
                  pool, iconv, ires, iskip, isk, sk, xx.in, ida, da, typeu, mix,
                  n11, n33, pa1, ikr, jkr, ikkk, ikkkm, cexp, interval, eps)

            if (z1$ier == 0) {
              iaddflg <- 1
              idx <- idx + 1
              caic <- c(caic, z1$aic)
              cexvar[[idx]] <- z1$exvar
              nrange[[idx]] <- z1$nrange
              cinterval[[idx]] <- z1$range
              ctbl[[idx]] <- z1$n
              ctbl.p[[idx]] <- z1$p
            } # if (z1$ier == 0) end

          } # if (iaddflg == 0) end 
        }  # if (i in 2:nsub) end
      }  # if (nsub > 2) end


      catdap2.out <- list(title = title, accuracy = xx, ires = ires,
           print.level = print.level, tway.table = cross.table, total = total,
           interval = interval, base.aic = baic, aic = aic,
           aic.order = aic.order, nsub = nsub, subset = list(nv = nv, ncc = ncc,
           aic = aaic, exv = cexp, vname = vnames),
           ctable = list(aic = caic, exvar = cexvar, nrange = nrange,
                         range = cinterval, n = ctbl, p = ctbl.p), 
           missing = nmiss)

      if (plot != 0)
        plot.catdap2(catdap2.out, plot, gray.shade)

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
      catdap2.out <- list(ier = ier)
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
    tmp <- as.factor(data)
    cdata <- as.numeric(tmp)
    pool <- 2
    xx <- 0.0
    dname <- levels(tmp)
  }

  return(list(pooling = pool, acc = xx, cdata = cdata, nmiss = nmiss,
              dname = dname, n1 = n1))

}

#================================================================

addaicm <- 
 function(iex2, nsamp, n, l, recode, iskip1, it, nov, item1, item2, pool,
 iconv, ires, iskip, isk, sk, xx, ida, da, typeu, nmiss, n11, n33, pa1, ikr,
 jkr, ikkk, ikkkm, cexp, interval, eps) {

#================================================================

    icl <- 1
    icl1 <- icl + 1
    nex <- length(cexp[[iex2]])
    icls <- c(nex+1, ires)
    for (i in 1:nex)
      icls <- c(icls, cexp[[iex2]][i])

    z <- .Fortran(C_catdap2m,
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
                  iab = integer(n*n11*n33),
                  totalc = integer(n11),
                  ttrr = integer(n*n33),
                  ab = double(n*n33),
                  iaa = integer(n*n33),
                  pa = double(n*n11*n33),
                  idata = integer(n),
                  ite = integer(n),
                  dx = double(n),
                  aaam = double(ikr),
                  caa = integer(ikr*jkr),
                  icaa = integer(ikr),
                  nnia = integer(ikr),
                  lk77 = integer(1),
                  morder = integer(ikr),
                  ibc = integer(n11*ikkk*icl1),
                  pbc = double(n11*ikkk*icl1),
                  aic1 = double(icl1),
                  iabse = integer(nov*n33*icl1),
                  baic = double(1),
                  as.integer(n11),
                  as.integer(n33),
                  as.integer(ikr),
                  as.integer(jkr),
                  as.integer(ikkk),
                  as.double(eps),
                  nrange = integer(nov*icl1),
                  ier = integer(3))

    ier <- z$ier
 
    if (ier[1] == 0 ) {
      iaa <- array(z$iaa, dim=c(n, n33))
      ite <- z$ite
      ibc <- array(z$ibc, dim=c(n11, ikkk, icl1))
      pbc <- array(z$pbc, dim=c(n11, ikkk, icl1))
      aic1 <- z$aic1
      iabse <- array(z$iabse, dim=c(nov, n33, icl1))
      nr <- array(z$nrange, dim=c(nov, icl1))
      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      i1 <- 1
      i2 <- ite[ires]

      intval <- list()
      np1 <- 1
      cint <- NULL
      exvar <- NULL
      exrange <- NULL
      nvv <- icls[1] - 1
      for (i in 1:nvv) {
         exv <- icls[i+2]
         exnr <- nr[i, 2]
         intv <- NULL
         j1 <- 1
         j2 <- ite[exv]
         if (xx[exv] != 0.0) {
            p <- iabse[(i+1), j1:(j2+1), 2]
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
         exrange <- c(exrange, exnr)
      }

      cinterval <- cint
      cexvar <- exvar
      nrange <- exrange
      ctbl <- ibc[1:np1, i1:i2, 2]
      ctbl.p <- pbc[1:np1, i1:i2, 2]
      caic <- aic1[2]

      return(list(aic = caic, exvar = cexvar, nrange = nrange,
                  range = cinterval, n = ctbl, p = ctbl.p, ier = 0))


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
          stop(" Working area for second contingency table is too short. pa1 can no longer set a larger value.")
        if ( ier[2] != -999)
          stop(gettextf(" Working area for contingency table is too short, try pa1== %i.", pa1n), domain = NA)
      }  # if (ier[1] == 2048) end
      return(list(ier = ier[1]))
    }
}

#==================================

print.catdap2 <- function(x, ...) {

#==================================

    eps <- 1.0e-10   # an error tolerance in the difference of AIC's
    dl18 <- "------------------"
    dl10 <- "----------"
    bl15 <- "               "
    bl8 <- "        "

    if (is.null(x$ier[1])) {
      n <- length(x$title)
      res <- x$ires
      lk <- x$nsub
      ktt <- max(x$subset$nv)

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
      cat("\n<< List of single explanatory variables (arranged in ascending order of AIC) >>\n\n")
#----------------------------------------------------------------------------#
      cat(sprintf(" Response variable : %s \t(base AIC = %8.2f)\n",
          x$title[res], x$base.aic))
      cat(rep(dl10,8), sep = "")
	  cat("\n")
      cat(rep(bl15, 2), "Number of\n", sep = "")
      cat(bl8, "      Explanatory     categories", bl15, "Difference\n", sep = "")
      cat(bl8, "      variables       of exp. var.    A I C      of AIC",
          bl8, "Weight\n",  sep = "")
      cat(rep(dl10,8), sep = "")
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

        cat(bl15, "    ", sep="")
        for (ii in i1:i2)
          cat(sprintf("%15s   ", dname[[res]][ii]))
        cat("          Total\n")
        cat(sprintf("( %15s )\n", x$title[j3]))
        for (ii in i1:(i2+2))
          cat(dl18, sep = "")
        cat("\n")
        for (k in j1:j2) {
          ptr <- 0
          if (x$total[[j3]][k] != 0)
            ptr <- (x$total[[j3]][k] * 100.) / x$total[[j3]][k]
          cat(sprintf("%8i          ", k))
          for (ii in i1:i2)
            cat(sprintf("%8i ( %5.1f )", x$tway.table[[j3]]$n[k,ii],
                x$tway.table[[j3]]$p[k,ii]))
          cat(sprintf("%8i ( %5.1f )\n", x$total[[j3]][k], ptr))
        }  # for (k in j1:j2) end

        for (ii in i1:(i2+2))
          cat(dl18, sep = "")
        cat("\n     Total        ")
        for (ii in i1:i2)
          cat(sprintf("%8i ( %5.1f )", x$total[[res]][ii], ptc[ii]))
        cat(sprintf("%8i ( %5.1f )\n", nsamp, ptt))

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
            lk4 <- x$subset$nv[ij]
            if (lk4 == j) {
              aaaa <- x$subset$aic[ij]
              if (ijk == 0)
                aaic1 <- x$subset$aic[ij]
              daic <- aaaa - aaic1
              w <- exp(-1. / 2 * daic)
              ijk <- ijk + 1

              if (ijk == 1) {
                cat(sprintf("\n Number of explanatory variables = %i\n", j))
                cat(rep(dl10, 8), sep = "")
                cat("\n", rep(bl15, 2), "Number of\n", sep = "")
                cat(bl8,
                    "      Explanatory     categories", bl15, " Difference\n", sep = "")
                cat(bl8,
                    "      variables       of exp. var.    A I C      of AIC",
                    bl8, "Weight\n", sep = "")
                cat(rep(dl10, 8), sep = "")
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
                cat(sprintf("%20s", x$subset$vname[[ij]][1]))
                cat(sprintf("     %8i", x$subset$ncc[ij]))
                cat(sprintf("     %8.2f", x$subset$aic[ij]))
                cat(sprintf("     %8.2f", daic))
                cat(sprintf("     %8.2f\n", w))
              }

              if (lk4 == 2 || lk4 > 2)
                for (i in 2:lk4)
                  cat(sprintf("     %20s\n", x$subset$vname[[ij]][i]))

            }  # if (lk4 == j) end
            aaic2 <- aaaa
          }  # for (ij in 1:lk) end
        }  # for (j in 1:ktt) end
      }  # if (x$print.level == 0) end

#-------------------------------------------------------------------#
      cat("\n<< Summary of subsets of explanatory variables >>\n")
#-------------------------------------------------------------------#

      cat(sprintf("\n Response variable : %s\n", x$title[res]))
      cat(rep(dl10, 8), sep = "")
      cat("\n", rep(bl15, 2), "Number of\n", sep = "")
      cat(bl8, "      Explanatory     categories", bl15, " Difference\n", sep = "")
      cat(bl8,
          "      variables       of exp. var.    A I C      of AIC",
          bl8, "Weight\n", sep = "")
      cat(rep(dl10, 8), sep = "")
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
        lk4 <- x$subset$nv[ij]
        if (lk4 > 0) {
          aaaa <- x$subset$aic[ij]
          if (ijk == 0 )
            aaic1 <- x$subset$aic[ij]
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
          cat(sprintf("%20s", x$subset$vname[[ij]][1]))
          cat(sprintf("     %8i", x$subset$ncc[ij]))
          cat(sprintf("     %8.2f", x$subset$aic[ij]))
          cat(sprintf("     %8.2f", daic))
          cat(sprintf("     %8.2f\n", w))

          if (lk4 == 2 || lk4 > 2)
            for (i in 2:lk4) {
              cat(sprintf("     %20s\n", x$subset$vname[[ij]][i]))
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
#--
            nrank <- 1
#--
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

      nctbl <- length(x$ctable$aic)
      icl <- nctbl - 1
      if ((icl != 0) && (x$ctable$aic[nctbl] == x$ctable$aic[nctbl-1]))
        icl <- icl - 1

      for (k in 1:nctbl) {
        lk4 <- 0
        if (length(x$ctable$exvar) > (k-1))
          lk4 <- length(x$ctable$exvar[[k]])
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
            lkk <- x$ctable$exvar[[k]][i]
            cat(sprintf(" X(%i) : %s\n", (i+1), x$title[lkk]))
            etitle <- c(etitle, x$title[lkk])
            j2[i] <- length(x$ctable$range[[k]][[i]])
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
        if (lk4 == 1)
          cat("    ", sep = "")

        for (i in i1:i2)
          cat(sprintf("           %3i    ", i))
        cat(bl8, " Total\n")
        for (ii in i1:(i2+1))
          cat(dl18, sep = "")
        cat(rep("----", max(2,lk4)), "\n", sep = "")

        idf <- TRUE
        if (lk4> 0) {
          nr <- x$ctable$nrange[[k]]
          np <- nr[1]
          if (np > 1 && lk4 > 1)
            for (i in 2:lk4)
              np <- np * nr[i]
          cidx <- matrix(, nrow=np, ncol=lk4)
          ipos <- lk4
          cidx[1, ] <- rep(1, lk4)

          for (ii in 2:np) {
            iline <- cidx[ii-1, ]
            kflg <- FALSE
            for (jj in lk4:1)
              if (kflg == FALSE) {
                jdx <- iline[jj] + 1
                if ( jdx < nr[jj] || jdx == nr[jj]) {
                  iline[jj] <- iline[jj] + 1
                  cidx[ii, ] <- iline
                  kflg <- TRUE
                } else {
                  iline[jj] <- 1
                }
              }
          }

          for (i in 1:np) {
            for (j in 1:lk4)
                cat(sprintf(" %2i ", cidx[i,j]))
            if (lk4 == 1)
              cat("    ")
            for (j in i1:i2)
              cat(sprintf("%8i ( %5.1f )", x$ctable$n[[k]][i,j],
                  x$ctable$p[[k]][i,j]))
            ibct <- sum(x$ctable$n[[k]][i, i1:i2])
            pbct <- sum(x$ctable$p[[k]][i, i1:i2])
            cat(sprintf("%8i ( %5.1f )\n", ibct, pbct))
          }  # end for (i in 1:np)

          for (ii in i1:(i2+1))
            cat(dl18)
          cat(rep("----", max(2, lk4)), "\n", sep = "")

        }   # if (lk4 > 0) end

        cat("   Total")
        if (lk4 > 2)
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
            lkk <- x$ctable$exvar[[k]][i]

            nmiss <- x$missing[[lkk]]
            ntype <- length(nmiss)
            if (ntype==1 && nmiss[1]==0)
              ntype <- 0

            ix <- i + 1
            print.CtableNote(ix, x$title[[lkk]], x$accuracy[[lkk]],
                             x$ctable$range[[k]][[i]], dname[[lkk]], j1,
                             j2[i], ntype, nmiss)

          }   #  for (i in 1:lk4) end
        }   # if (lk4 > 0) end

        cat(sprintf("\nAIC = %8.2f\n", x$ctable$aic[k]))
        cat(sprintf("base AIC = %8.2f\n", x$base.aic))

      } # for (k in  1:ntbl) end

      nbest <- 1
      if (lk > 1) { 
        aicm1 <- x$subset$aic[1]
        for (i in 2:lk) {
          aicm2 <- x$subset$aic[i]
#    2020/09/06
          if (aicm1 == 0 && aicm2 == 0)
            break

          daic <- abs(aicm2 - aicm1) / max(abs(aicm2), abs(aicm1))
          if (daic > eps)
            break
          nbest <- nbest + 1
        }
        if (nbest == 3)
          cat("<NOTE> There is another subset with the minimum AIC.\n\n")
        if (nbest > 3)
          cat(sprintf("<NOTE> There are another %3i subsets with the minimum AIC.\n", nbest-2))
      }
    } else {
      ier <- x$ier[1]
      eval <- x$ier[2]

      if (ier == 2002) {
        stop(" Working area for multi-dimensional table is not enough. Try larger pa2.")
      } else if (ier == 2003) {
        stop(" Working area for multi-dimensional table is not enough. Try larger pa3.")
      } else if (ier == 2048) {
        if (eval != -999)
          stop(gettextf(" Working area for contingency table is too short, try pa1= %i.", eval), domain = NA)
        if (eval == -999)
          stop(" Working area for contingency table is too short. pa1 can no longer set a larger value.")
      } else if (ier == 2037) {
        stop(" the program catdap cannot deal with data sets where the number of non-zero frequency categories of the response variables is less than 2.")
      } else if (ier == 2035) {
        stop(" The value of variable is beyond the interval specified in 'min' and 'max'.")
      } else if (ier == 2588) {
        stop(" lk5 > n-1 ")
      } else if (ier == 650) {
        stop(" The value of 'nvar' is too small for the additional analysis.")
      }
    }

}


#============================================

plot.catdap2 <- function(x, plot.type, gray.shade, ...) {

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

    nctbl <- length(x$ctable$aic)
    icl <- nctbl - 1
    if ((icl != 0) && (x$ctable$aic[nctbl] == x$ctable$aic[nctbl-1]))
      icl <- icl - 1

    icflg <- 1
    for (k in 1:nctbl) {
      ctable.dname <- list()
      if (x$accuracy[[res]] == 0.0)
        ctable.dname[[1]] <- dname[[res]]
      if (x$accuracy[[res]] != 0.0)
        ctable.dname[[1]] <- c(i1:i2)

      lk4 <- 0
      if (length(x$ctable$exvar) > (k-1))
        lk4 <- length(x$ctable$exvar[[k]])
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
          lkk <- x$ctable$exvar[[k]][i]
          etitle <- c(etitle, x$title[lkk])
          j2[i] <- length(x$ctable$range[[k]][[i]])
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

      plot.mosaic(plot.type, gray.shade, res, icl, k, mp, idm, x$title, etitle,
                  dname, x$ctable$n, ctable.dname, nc, x$aic.order, x$aic,
                  x$tway.table, x$total, icflg)

    } # for (k in  1:ntbl) end


}


#==============================================================================

print.Note <- function(title, accuracy, interval, dname, i1, i2, ntype, nmiss) {

#==============================================================================

      cat("\n<Note>\n")
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
    cat("\n")
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
