#========================================================================

catdap1 <- function(cdata, response.names = NULL, plot = 1, gray.shade = FALSE,
                    ask = TRUE) {

#========================================================================

    if (plot != 1 && plot != 2)
      plot <- 0

#   cdata                    # original data
    n <- dim(cdata)[2]       # total number of variable
    nsamp <- dim(cdata)[1]   # sample size

#   title                   # variable names
    title <- names(cdata)
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

    catdata <- t(cdata)
    for (i in 1:n)
      if (is.numeric(catdata[i, ]) == FALSE && is.factor(catdata[i, ]) == FALSE)
        catdata[i, ] <- as.factor(catdata[i, ])


    dname <- list()
    for (i in 1:n)
      if (is.numeric(cdata[, i]) == FALSE) 
        dname[[i]] <- levels(ordered(cdata[, i]))

    catdap1.out <- catdap01(catdata, title, dname, response.names) 

    if (plot != 0)
      plot.catdap1(catdap1.out, plot, gray.shade, ask)

    class(catdap1.out) <- "catdap1"
    return(catdap1.out)
}

#==========================================================================

catdap1c <- function(ctable, response.names = NULL, plot = 1, gray.shade = FALSE,
                     ask = TRUE) {

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
      plot.catdap1(catdap1c.out, plot, gray.shade, ask)

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
        stop(" Response variable name is wrong.")
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
    for (i in 1:n) {
      item1[i] <- as.integer(min(cdata[i, ]))
      item2[i] <- as.integer(max(cdata[i, ]))
    }
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
    
    z <- .Fortran(C_catdap1,
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
                  as.integer(cdata),
                  nc = integer(n),
                  ia = integer(n*l*n4*n4),
                  p = double(n*l*n4*n4),
                  total = integer(n*n4),
                  aic = double(n*l),
                  ord = integer((n-1)*l))

    nc <- z$nc
    ia <- array(z$ia, dim=c(n4, n4, l, n))
    p <- array(z$p, dim=c(n4, n4, l, n))
    total <- array(z$total, dim=c(n4, n))
    aic <- array(z$aic, dim=c(l, n))
    for (i in 1:l)
      aic[i, ires[i]] <- NA
    ord <- array(z$ord, dim=c(l, n-1))
    ia <- aperm(ia, c(2, 1, 3, 4))
    p <- aperm(p, c(2, 1, 3, 4))

    cross.table <- list()
    for (i in 1:l) {
      i1 <- 1
      i2 <- nc[ires[i]]
      cross.table[[i]] <- list()
      for (j in 1:n) {
        cross.table[[i]][[j]] <- list()
        if (j != ires[i]){
          cross.table[[i]][[j]]$res <- title[ires[i]]
          cross.table[[i]][[j]]$exv <- title[j]
          j1 <- 1
          j2 <- nc[j]
          cross.table[[i]][[j]]$n <- array(ia[j1:j2, i1:i2, i, j], dim=c(j2,i2),
            dimnames = list(dname[[j]], dname[[ires[i]]]))
          cross.table[[i]][[j]]$p <- array(p[j1:j2, i1:i2, i, j], dim=c(j2,i2))
        } else {
          cross.table[[i]][[j]]$res <- title[ires[i]]
          cross.table[[i]][[j]]$exv <- NULL
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
    dl22 <- "----------------------"
    dl10 <- "----------"
    bl14 <- "              "
    bl8 <- "        "

#--------------------------------------------------------------#
    cat("\n < Summary of AIC's for the two-way tables >\n\n")
#--------------------------------------------------------------#
    cat("         Explanatory\n")
    cat("           variables")
    for (i in 1:l)
      cat(sprintf("  %20s", x$tway.table[[i]][[1]]$res))
    cat("\n", rep(dl10, 2), rep(dl22, l), "\n", sep="")
    for (jj in 1:n) {
      cat(sprintf("%20s", x$title[jj]))
      for (j in 1:l) {
        if (is.na(x$aic[j,jj]) == TRUE) {
          cat(bl14, "    -  ")
          ires[j] <- jj
        } else {
          if (is.null(x$aic[j,jj]) == FALSE)
            cat(bl14, sprintf("%8.2f", x$aic[j,jj]), sep = "")
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
      cat("\n\n < List of explanatory variables arranged in ascending order of AIC >\n")
#------------------------------------------------------------------------------#
      cat(sprintf("\n Response variable  : %s\n\n", x$title[i3]))
      cat("  No.         Explanatory     Number of       A I C     Difference\n")
      cat(bl14, "variables       categories", bl8, bl8, "of AIC\n", sep = "")
      cat(rep(dl10,7), sep = "")
      cat("\n")

      as <- 0.0
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        aic <- x$aic[i, nv]
        if (j != 1)
          as <- aic - as
        cat(sprintf("%5i%20s     %8i     %8.2f     %8.2f\n", j,x$title[nv], item[nv],
            aic, as))
        as <- aic
      }  # for (j in 1:nc) end
#--------------------------------------------------------------------------#
      cat("\n\n < Two-way tables arranged in ascending order of AIC >\n\n")
#--------------------------------------------------------------------------#

      cat(sprintf("\t\t\t\t\t( %s )\n", x$title[i3]))
      nsamp <- 0
      ptt <- 0.0
      pt <- rep(0, i2)
      cat(bl8)
      for (ii in i1:i2) {
        cat(sprintf("%18s", dname[[i3]][ii]))
        nsamp <- nsamp + x$total[[i3]][ii]
      }
      for (ii in i1:i2) {
        pt[ii] <- x$total[[i3]][ii] / nsamp * 100.0
        ptt <- ptt + pt[ii]
      }

      cat("             Total\n")
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        j1 <- 1
        j2 <- item[nv]
        cat(sprintf("( %s )\n", x$title[nv]))
        for( k in j1:j2) {
          cat(sprintf("%10s", dname[[nv]][k]))
          tc <- 0
          for (ii in i1:i2)
            cat(sprintf("% 8i ( %5.1f )", x$tway.table[[i]][[nv]]$n[k,ii],
                x$tway.table[[i]][[nv]]$p[k,ii]))
          cat(sprintf("% 8i ( %5.1f )\n", x$total[[nv]][k], ptt))
        }
        cat("     Total")
        for (ii in i1:i2)
          cat(sprintf("%8i ( %5.1f )", x$total[[i3]][ii], pt[ii]))
        cat(sprintf("%8i ( %5.1f )\n\n", nsamp, ptt))
      }  # for (j in 1:nc) end
    }  # for (i in 1:l) end

    if (x$ier == 3001)
      warning(" some response variable name is wrong")

}


#=================================================

plot.catdap1 <- function(x, plot.type, gray.shade, ask, ...) {

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
        plot.single1(dname, ires[i], item, x$aic.order[i, ], x$aic[i, ],
                     x$tway.table[[i]], x$total, gray.shade, old.par, par.ask)
      if (plot.type == 2)
        plot.single2(dname, ires[i], item, x$aic.order[i, ], x$aic[i, ],
                     x$tway.table[[i]], x$total, gray.shade, old.par, par.ask)
      par.ask <- ask
    } 

    i1 <- 1
    i3 <- ires[1]
    i2 <- item[i3]
    nsamp <- 0
    for (ii in i1:i2)
      nsamp <- nsamp + x$total[[i3]][ii]
    plot.shade(x$title, x$aic, ires, nsamp, gray.shade, old.par, ask)
}

#=================================================

plot.shade <- function(title, aic, ires, nsamp, gray.shade, old.par, ask) {

#=================================================

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

  mtitle1 <- "The shaded style display of all the AIC's\n\nColor level : "
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

  pcol <- brewer.pal(11, "RdGy")
  if (gray.shade == TRUE) {
    scol <- pcol[c(11,10,9,8,6)]
  } else {
    scol <- pcol[2:6]
  }

  image(x, y, aics2, col = scol, breaks = bp, axes = FALSE, xlab = "", ylab = "",
        main = main.title, sub = "explanatory variable", cex.main = 0.8,
        cex.sub = 0.8)

  legend("topright", c("5","4","3","2","1"), pch = 15, col = scol, bg = "snow2",
         cex = 0.8, pt.cex = 1.1)

  axis(1, 1:d2, explanatory, las = 3)
  axis(2, 1:d1, response, las = 1)
  box()

  par(old.par)
}


