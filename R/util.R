#======================================================================

# plot function (plot.single1, plot.single2, plot.mosaic)

#======================================================================

plot.single1 <- function(dname, res, item, aic.order, aic, tway.table,
                         total, gray.shade, old.par, ask) {

    n <- length(aic)
    i2 <- item[res]
    tt <- tway.table

    max.nc <- 4
    if (n < max.nc + 1) {
      nc <- n
      nr <- 1
    } else if (n < (max.nc * 2 + 1)) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- max.nc
      nr <- 2
    }
    m <- nc * nr
    insetd <- -0.4
    if (nc > 3) insetd <- -0.6

    new.mar <- c(4,3,2,3)    # default : c(5,4,4,2)
    new.oma <- c(3,3,2,3)    # default : c(3,3,3,3)

    if (gray.shade == T) {
      col <- TRUE
    } else {
      pcol <- brewer.pal(12, "Paired")
      col <- c(pcol[c(1,3,5,7,9,11)], pcol[c(2,4,6,8,10,12)])
    }

    nplot <- 1
    irflag <- 0
    nw <- 0
    if (ask ==TRUE)
      par(ask=TRUE)

    j1 <- 1
 	ngraph <- 0

    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tt[[iex]]$n
      resname <- tt[[iex]]$res
      exname <- tt[[iex]]$exv

      if (aic[iex] > 0 && irflag == 0) {
        if (nplot%%m == 1) {
          if (j > 1 || nw > 0)
            if (.Device != "null device" || ask == FALSE)
              dev.new()
          par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
          nw <- nw + 1
        }
        irflag <- 1
        y <- total[[res]]
        y <- as.array(y)
        xlabel <- paste(resname, " Total ", total[[res]][1])
        if (i2 > 1)
          for (ii in 2:i2)
            xlabel <- paste(xlabel, ":", total[[res]][ii])
        dimnames(y) <- list(dname[[res]])
        mosaicplot(y, color = col, main = "", xlab = "", ylab = "AIC= 0.0\n ",
                   off = 0)
        title(xlab = xlabel, line = 0.5)
				   
        nplot <- nplot + 1
        ngraph <- ngraph + 1

        if ((ngraph %% nc) == 0)
          legend("bottom", legend = dname[[res]], col = col, pt.cex = 1.5,
                 pch = 15, bty = "n", horiz = TRUE, inset = insetd)

      }

      if (nplot%%m == 1) {
        if (j > 1 || nw > 0)
          if (.Device != "null device" || ask == FALSE)
            dev.new()
        par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
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
      mosaicplot(y, color = col, main = "", xlab = "",
                 ylab = paste("AIC=", aaic, "\n", exname),
                 dir = c("h", "v"), off = c(10, 0))
      title(xlab = resname, line = 0.5)
				 
      nplot <- nplot + 1
      ngraph <- ngraph + 1

      if ((ngraph %% nc) == 0 || ngraph == n)
        legend("bottom", legend = dname[[res]], col = col, pt.cex = 1.5,
               pch = 15, bty = "n", horiz = TRUE, inset = insetd)

      if (nw == 1)
        mtext(text = "Single Explanatory Models in ascending order of AIC",
              side = 3, outer = TRUE, cex = 1.2)

    }  # for (j in 1:(n-1)) end

    if (irflag == 0) {
      if (nplot%%m == 1) {
        if (j > 1 || nw > 0)
          if (.Device != "null device" || ask == FALSE)
            dev.new()
        par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
      }

      y <- total[[res]]
      y <- as.array(y)
      xlabel <- paste(resname, " Total ", total[[res]][1])
      if (i2 > 1)
        for (ii in 2:i2)
          xlabel <- paste(xlabel, ":", total[[res]][ii])
      dimnames(y) <- list(dname[[res]])
      mosaicplot(y, color = col, main = "", xlab = "", ylab = "AIC= 0.0\n ",
                 off = 0)
      title(xlab = xlabel, line = 0.5)

      nplot = nplot + 1
      ngraph <- ngraph + 1
      if ((ngraph %% nc) == 0 || ngraph == n)
        legend("bottom", legend = dname[[res]], col = col, pt.cex = 1.5,
               pch = 15, bty = "n", horiz = TRUE, inset = insetd)
    }

    par(old.par)

}

plot.single2 <- function(dname, res, item, aic.order, aic, tway.table,
                         total, gray.shade, old.par, ask) {

    n <- length(aic)
    i2 <- item[res]
    tt <- tway.table

    max.nc <- 4
    if (n < max.nc + 1) {
      nc <- n
      nr <- 1
    } else if (n < (max.nc * 2 + 1)) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- max.nc
      nr <- 2
    }
    m <- nc * nr

    new.mar <- c(4,3,2,3)    # default : c(5,4,4,2)
    new.oma <- c(3,3,2,3)    # default : c(3,3,3,3)

    if (gray.shade == T) {
      col <- TRUE
    } else {
      pcol <- brewer.pal(12, "Paired")
      col <- c(pcol[c(1,3,5,7,9,11)], pcol[c(2,4,6,8,10,12)])
    }

    nplot <- 1
    irflag <- 0
    nw <- 0
    if (ask ==TRUE)
      par(ask=TRUE)

    j1 <- 1
    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tt[[iex]]$n
      resname <- tt[[iex]]$res
      exname <- tt[[iex]]$exv

      if (aic[iex] > 0 && irflag == 0) {
        if (nplot%%m == 1) {
          if (j > 1 || nw > 0)
            if (.Device != "null device" || ask == FALSE)
              dev.new()
          par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
          nw <- nw + 1
        }
        irflag <- 1
        yy <- total[[res]]
        y <- matrix(yy, nrow = length(yy))
        xlabel <- paste(resname, " Total ", total[[res]][1])
        if (i2 > 1)
          for (ii in 2:i2)
            xlabel <- paste(xlabel, ":", total[[res]][ii])
        dimnames(y) <- list(dname[[res]])
        colnames(y) <- rep("")
        mosaicplot(y, color = col, main = "", xlab = "", ylab = "AIC= 0.0\n ")
        title(xlab = xlabel, line = 0.5)

        nplot <- nplot + 1

      }

      if (nplot %% m == 1) {
        if (j > 1 || nw > 0)
        if (.Device != "null device" || ask == FALSE)
          dev.new()
        par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
        nw <- nw + 1
      }

      j2 <- item[iex]
      y <- array(0, dim=c(i2, j2))
      for (ii in 1:i2) 
        for (k in j1:j2) y[ii, k] <- ctbl[k, ii]
      dimnames(y)[1] <- list(dimnames(ctbl)[[2]])
      dimnames(y)[2] <- list(dimnames(ctbl)[[1]])
      aaic <- round(aic[iex], 2)
      mosaicplot(y, color = col, main = "", xlab = "",
                 ylab = paste("AIC=", aaic, "\n", exname))
      title(xlab = paste(resname), line = 0.5)

      nplot <- nplot + 1
      legend(x = 1, y = 0.5, legend = dimnames(ctbl)[[1]], col = col,
             pt.cex = 1.5, pch = 15, bty = "n")

      if (nw == 1)
        mtext(text = "Single Explanatory Models in ascending order of AIC",
              side = 3, outer = TRUE, cex = 1.2)

    } 

    if (irflag == 0) {
      if (nplot%%m == 1) {
        if (j > 1 || nw > 0)
          if (.Device != "null device" || ask == FALSE)
            dev.new()
        par(mfcol=c(nc, nr), mar = new.mar, oma = new.oma, xpd = T)
      }

      yy <- total[[res]]
      y <- matrix(yy, nrow = length(yy))
      xlabel <- paste(resname, " Total ", total[[res]][1])
      if (i2 > 1)
        for (ii in 2:i2)
          xlabel <- paste(xlabel, ":", total[[res]][ii])
      dimnames(y) <- list(dname[[res]])
      colnames(y) <- rep("")
      mosaicplot(y, color = col, main = "", xlab = "", ylab = "AIC= 0.0\n ")
      title(xlab = xlabel, line = 0.5)
   }

    par(old.par)

}

plot.mosaic <- function(iplt, gray.shade, res, icl, k, mp, idm, title, etitle,
                        dname, ctable, ctable.dname, nc, aic.order, aic,
                        tway.table, total, iflg)

{
  old.par <- par(no.readonly = TRUE)
  ask <- FALSE

  if (k == 1) {
    if (iplt == 1)
      plot.single1(dname, res, nc, aic.order, aic, tway.table, total,
                   gray.shade, old.par, ask)
    if (iplt == 2)
      plot.single2(dname, res, nc, aic.order, aic, tway.table, total,
                   gray.shade, old.par, ask)
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
    if (gray.shade == TRUE ) {
      mosaicplot(t, color = TRUE, main = mtitle, xlab = paste(xtitle),
                 ylab = paste(ytitle), dir = split, off = offset)
    } else {
      pcol <- brewer.pal(12, "Paired")
      col <- c(pcol[c(1,3,5,7,9,11)], pcol[c(2,4,6,8,10,12)])
      mosaicplot(t, color = col, main = mtitle, xlab = paste(xtitle),
                 ylab = paste(ytitle), dir = split, off = offset)
    }

    if (iplt == 1) {
      par(xpd = T)
 	  dname.res <- ctable.dname[[1]]
      legend("bottom", legend = dname.res, col = col, pt.cex = 1.5, pch = 15,
             bty = "n", horiz = TRUE, inset = -0.15)
    }

  } # if (ifg > 0) end

    par(old.par)
}


#=====================================================================

Barplot2WayTable <- function (x, exvar = NULL, gray.shade = FALSE) {

#=====================================================================

  if (is(x, "catdap1") == FALSE && is(x, "catdap2") == FALSE)
    stop(" Object class is invalid.")

  tt <- x$tway.table
  vname <- x$title
  old.par <- par(no.readonly = TRUE)

  if (gray.shade == FALSE) {
    pcol2 <- brewer.pal(8, "Set2")
    pcol3 <- brewer.pal(12, "Set3")
    gcol <- c(pcol2, pcol3)
  }

  if (is(x, "catdap1")) {  # catdap1() or catdap1c()

    nres <- length(tt)
    nv <- length(x$title)
    if (is.null(exvar) == TRUE)
      exvar <- x$title
    nexv <- length(exvar)
    n <- min(nexv, nv-1)

    max.nc <- 4
    if (n < max.nc + 1) {
      nc <- n
      nr <- 1
    } else if (n < (max.nc * 2 + 1)) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- max.nc
      nr <- 2
    }
    m <- nc * nr

    par(mfcol = c(nc, nr), mar = c(5,4,3,5)+0.1)

    for (i in 1:nres) {

      if (i != 1)
        dev.new()
      par(mfcol = c(nc, nr), mar = c(5,4,3,5)+0.1)
      nplot <- 0

      for (j in 1:nexv) {
        ex <- 0
        for (k in 1:nv)
          if(exvar[j] == vname[k]) {
            ex <- k
            break
        }
        if (ex == 0)
          stop(gettextf(" %s is invalid for an explanatory variable name.", exvar[j]), domain = NA)

        if (is.null(tt[[i]][[ex]]$exv) == FALSE) {
          if (nplot != 0 && (nplot+1)%%m == 1) {
             dev.new()
             par(mfcol = c(nc, nr), mar = c(5,4,3,5)+0.1)
          }
          resname <- tt[[i]][[ex]]$res
          exname <- tt[[i]][[ex]]$exv
          h <- tt[[i]][[ex]]$n

          nr <- length(dimnames(h)[[1]])
          if (gray.shade == TRUE) {
            ccol <- seq(2, 10, by=9/(nr+1))
            gcol <- c(gray(ccol*0.1))
            gcol <- gcol[1:nr]
          }

          barplot(h, col = gcol, names.arg = dimnames(h)[[2]], legend = dimnames(h)[[1]],
                  args.legend = list(x="topright", title=exname, cex=0.8, inset=c(-0.1,-0.1)),
                  space = 0, xlab = resname)
          nplot <- nplot + 1
        }
      }  # for (j in 1:nexv) end
    }  # for (i in 1:nres) end

  } else if (is(x, "catdap2")) {  # catdap2()

    nv <- length(x$title)
    if (is.null(exvar) == TRUE)
      exvar <- x$title
    nexv <- length(exvar)
    n <- min(nexv, nv-1)
    interval <- x$interval

    max.nc <- 4
    if (n < max.nc + 1) {
      nc <- n
      nr <- 1
    } else if (n < (max.nc * 2 + 1)) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- max.nc
      nr <- 2
    }
    m <- nc * nr
    par(mfcol = c(nc, nr), mar = c(5,4,3,5)+0.1)

    nplot <- 0
    for (j in 1:nexv) {
      ex <- 0
      for (i in 1:nv)
        if(vname[i] == exvar[j]) {
          ex <- i
          break
        }
      if (ex == 0)
        stop(gettextf(" %s is invalid for an explanatory variable name.", exvar[j]), domain = NA)

      if (is.na(tt[[ex]][1]) == FALSE) {
        resname <- tt[[ex]]$res
        res <- 0
        for (i in 1:nv)
          if(vname[i] == resname) {
            res <- i
            break
          }

        if (nplot != 0 && (nplot+1)%%m == 1) {
          dev.new()
          par(mfcol = c(nc, nr), mar = c(5,4,3,5)+0.1)
        }

        h <- tt[[ex]]$n
        nr <- length(dimnames(h)[[1]])

        if (gray.shade == TRUE) {
          ccol <- seq(2, 10, by=9/(nr+1))
          gcol <- c(gray(ccol*0.1))
          gcol <- gcol[1:nr]
        }

        if (j == 1) {
          nrres <- length(dimnames(h)[[2]])
          resint <- dimnames(h)[[2]]
          resname <- tt[[ex]]$res
        }

        nint <- length(interval[[ex]])
        nmiss <- 0
        if (length(x$missing[[ex]]) > 0 && x$missing[[ex]][1] != 0)
          nmiss <- length(x$missing[[ex]])

        if (nr == nint) {
          barplot(h, col = gcol, names.arg = resint, legend = dimnames(h)[[1]],
                  args.legend = list(x="topright", title=vname[ex], cex=0.8, inset=c(-0.1,-0.1)),
                  space = 0, xlab = resname)
        } else {
          itext <- NULL
          for (k in 1:(nr-nmiss))
            itext[[k]] <- paste(interval[[ex]][k], "-", interval[[ex]][k+1])
          if (nmiss > 0)
            for (k in 1:nmiss)
              itext[[nr-nmiss+k]] <- paste("missing-", x$missing[[ex]][k])

          barplot(h, col = gcol, names.arg = resint, legend = itext,
                  args.legend = list(x="topright", title=vname[ex], cex=0.8, inset=c(-0.1,-0.1)),
                  space = 0, xlab = resname)
        }
        nplot <- nplot + 1
      }
    }  # for (j in 1:nexv) end

    nint <- length(interval[[res]])
    if (nint > nrres) {
      cat("\n<Note>\n")
      cat(sprintf(" %s\n", resname))
      cat("\tcategory    \tvalue range\n")
      for (k in 1:nrres)
        cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", k,
                interval[[res]][k], interval[[res]][k+1]))
    }
  }
  par(old.par)
}
