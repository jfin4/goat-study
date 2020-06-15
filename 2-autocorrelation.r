# setup ----------------------------------------------------------------
rm(list = ls())
library(ape)
set.seed(0)

# import dat ----------------------------------------------------------
dat <- read.csv("data-long.csv")
names(dat) <- c("can", "bout", "plot", "tmnt", "tsct", "cat", "occ")
dat$can <- ifelse(dat$can == "barloy", 2, 1)
dat$tmnt <- ifelse(dat$tmnt == "control", 0, 1)
dat$tid <- # give uniqe id to each transect-bout
    with(dat, paste0(can, bout, plot, tmnt, tsct))
dat$pt <- rep(1:40, len = nrow(dat))


# define functions --------------------------------------------------
get_i <- function(tran, mlag) {
    # get Moran's I of a transect at a given lag distance
    Moran.I(tran, mlag)$observed
}
get_is <- function(lag, bdata, mdist) {
    # get Moran's Is at each lag distance of each (resampled) transect-bout
    # assumes data comes in "box" format, which is a stack of tables. 
    #   data are occurrences at points along transects (0, 1)
    #   the bottom stack is observed data, 
    #       rows are transect-bouts,   
    #       columns are points along the transects-bouts
    #   the upper stacks are resampled transect-bout data
    mlag <- # create lag-weighted distance matrix
        ifelse(mdist == lag, 1, 0)
    # returns a table. Has same format as observed data table, except
    #   columns are lag distances instead of transect points
    apply(bdata, c(1, 3), get_i, mlag)
}


# number of iterations -------------------------------------------------
# 1000 takes a long time to run. Load image below to skip computation.
nits <- 1000

# permutate Is ---------------------------------------------------------
cats <- c("bapi", "bu", "ag", "onh", "onn", "inv")
lboxis <- lapply(1:length(cats), function(x) NA)
for (i in 1:length(cats)) {
    # i <- 1
    cat <- cats[i]
    datcat <- dat[dat$cat == cat, ]
    # create matrix or transect-bout data
    #     rows are transect-bouts
    #     columns are points along transects
    trans <- tapply(datcat$occ, list(datcat$tid, datcat$pt), function(x) x)
    trans <- # remove transect-bouts with 0 variance or Moran.I() will choke
        trans[apply(trans, 1, var) != 0, ]
    nrows <- nrow(trans)
    ncols <- ncol(trans)
    nlags <- # weirdness happens when nlags >= half the no. transect points
        ((ncol(trans) / 2) - 1 )
    mdist <- as.matrix(dist(1:ncols))
    # create stack of matrices of resampled transect-bout data
    #     rows are transect-bouts
    #     columns are points along transects
    #     layers are iterations (n = nits)
    reps <- replicate(nits, t(apply(trans, 1, sample, ncols)))
    boxcat <- # create the "box"
        array(NA, dim = c(nrows, ncols, nits + 1))
    boxcat[, , 1] <- # bottom layer is observed data
        trans
    boxcat[, , -1] <- # upper layers are resampled data
        reps
    # compute Moran's I for each (resampled) transect-bout at each lag distance
    # returns a list of tables
    lis <- lapply(0:(nlags + 1), get_is, bdata = boxcat, mdist = mdist)
    # Stack each element of the list into another box
    boxis <- array(unlist(lis), dim = c(nrows, nits + 1, nlags))
    # permutate box to same format where 
    #     rows are transect-bouts 
    #     columns are lag distances 
    #     layers are observed/resampled
    lboxis[[i]] <- aperm(boxis, c(1, 3, 2))
}


# calculate ranks
lranks <- lapply(1:length(cats), function(x) NA)
for (i in 1:length(cats)) {
    boxis <- lboxis[[i]]
    # "collapse" box into table of proportions where
    #   a datum is the proportion of resampled data "above" an observed datum 
    #     that is greater than or equal to that observed datum
    #   rows are transect-bouts
    #   columns are lag distances
    lranks[[i]] <- apply(boxis, c(1, 2), function(k) sum(k[1] >= k[-1]) / nits)
}

# calculate median ranks
lmedians <- lapply(1:length(cats), function(x) NA)
for (i in 1:length(cats)) {
    ranks <- lranks[[i]]
    # "collapse" table into vector of medias 
    lmedians[[i]] <- apply(ranks, 2, median)
}

# save r data 
# save.image("./autocorrelation.rdata")

# load r data #
load("./autocorrelation.rdata")


# plot
catnams <- c("Coyote Brush", "Bunch Grass", "Annual Grass",
             "Other Native Herb", "Other Non-native", "Invasive")
colors <- c("#996600", "#669900", "#ffcc00", 
          "#99cc00", "#ff9900", "#ff0000")
x <- 0.2
xsleg <- c(rep(x, 6))
ysleg <- rev(seq(0.03, 0.27, length.out = 6))
x1sleg <- rep(x + 5.1, 6)
x2sleg <- rep(x + 6.7, 6)
lwd <- 2
ltys <- c("solid", "a2", "8212", "621212", "42121212", "12")

# plot 1
data <- lboxis[[1]]
transect <- data[1, , 1]
treps <-  t(data[1, , -1])
png('./figures/correlation-1.png',
    width = 6, height = 4, units = "in", res = 300)
par(las = 1, mar = c(4, 5, 2, 1), 
    cex = 0.8, cex.axis = 0.9,
    xaxs = "i", yaxs = "i")
plot(c(1, 20), c(-1, 1),  type="n", axes = F, xlab = "",
     xlim = c(1, nlags + 2), ylim = c(-1, 1), 
     ylab = "Moran's I")
title(xlab = "Lag Distance (m)", line =2)
box(bty = "l")
boxplot(treps, add = T, border = "gray", col = "gray", outline = F, 
        outwex = 2.5, staplewex = 0, whisklty = "solid", range = 0, 
        boxwex = 0.7, axes = F)
lines(1:nlags, transect, lwd = 1, col = colors[1])
axis(1, at = c(1, 5, 9, 13, 17, 21), labels = c(0, 1:5))
axis(2)
lines(c(1.5, 2.3), c(-0.6, -0.6), lwd = 1, col = colors[1])
t1 <- expression(paste(italic(I), 
                      " of example transect at lag distance"))
text(2.6, -0.6, t1, adj = c(0, 0.5))
polygon(c(1.5, 1.5, 2.3, 2.3), c(-0.85, -0.75, -0.75, -0.85), 
        col = "gray", border = "white")
lines(c(1.9, 1.9), c(-0.7, -0.9), col = "gray")
t2 <- expression(paste("IQR and range of ", 
                      italic(I), 
                      " of resampled transect (n = 1000) at lag distance"))
text(2.6, -0.8, t2, adj = c(0, 0.5))
dev.off()

# plot 2
ranks <- lranks[[1]]
rank25 <- apply(ranks, 2, quantile, probs = 0.25)
rank75 <- apply(ranks, 2, quantile, probs = 0.75)
png(paste0('./figures/correlation-2.png'),
    width = 6, height = 4, units = "in", res = 300)
par(las = 1, mar = c(4, 5, 2, 1), 
    cex = 0.8, cex.axis = 0.9,
    xaxs = "i", yaxs = "i")
plot(lmedians[[1]], type = "n", axes = F, xlab = "",
     xlim = c(0, nlags + 1), ylim = c(0, 1), ylab = "")
title(xlab = "Lag Distance (m)", line = 2)
title(ylab = "Median Rank Relative to Resampled", line = 3)
lines(0:(nlags + 1), rep(0.5, nlags + 2), lwd = 0.5)
ltys <- c("solid", "a2", "8212", "621212", "42121212", "12")
i <- 1
medians <- lmedians[[i]]
lines(0:(nlags - 1), medians, lwd = lwd, col = colors[i], lty = ltys[i])
lines(0:(nlags - 1), rank25, lty = 3, col = colors[i])
lines(0:(nlags - 1), rank75, lty = 3, col = colors[i])
axis(1, at = c(0, 4, 8, 12, 16, nlags + 1), labels = c(0, 1:4, 5))
axis(2, at = c(0, 0.25, 0.5, 0.75, 1 ))
lines(c(0.5, 1.3), c(0.10, 0.10), lwd = 2, col = colors[1])
t1 <- expression(paste("Median Rank of Observed ",
                       italic(I), 
                      " Relative to Resampled"))
text(1.6, 0.10, t1, adj = c(0, 0.5))
lines(c(0.5, 1.3), c(0.05, 0.05), lty = 3, col = colors[1])
t2 <- expression(paste("25th/75th Percentile Rank of Observed ",
                       italic(I), "..."))
text(1.6, 0.05, t2, adj = c(0, 0.5))
dev.off()

# plot 3 -----------------------------------------------------------------
png(paste0('./figures/correlation-3.png'),
    width = 6, height = 4, units = "in", res = 300)
par(las = 1, mar = c(4, 5, 2, 1), 
    xaxs = "i", yaxs = "i",
    cex = 0.8, cex.axis = 0.9)
plot(lmedians[[1]], type = "n", axes = F, xlab = "",
     xlim = c(0, nlags + 1), ylim = c(0, 1), 
     ylab = "Median Rank Relative to Resampled")
title(xlab = "Lag Distance (m)", line = 2)
lines(0:(nlags + 1), rep(0.5, nlags + 2), lwd = 0.5)
for (i in length(cats):1) {
    medians <- lmedians[[i]]
    lines(0:(nlags - 1), medians, lwd = lwd, col = colors[i], lty = ltys[i])
    axis(1, at = c(0, 4, 8, 12, 16, nlags + 1), labels = c(0, 1:4, 5))
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1 ))
    text(x = xsleg[i], y = ysleg[i], labels = catnams[i], adj = c(0, 0.5))
    lines(x = c(x1sleg[i], x2sleg[i]), y = c(ysleg[i], ysleg[i]), 
          lwd = lwd, col = colors[i], lty = ltys[i])
}
dev.off()
