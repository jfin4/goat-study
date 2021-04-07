# setup ------------------------------------------------------------------------
rm(list = ls())
library(blme)

# import/prepare data ----------------------------------------------------------
data <- read.csv('./data-long.csv',
                 colClasses = "character")
names(data)[match(c("canyon", "treatment", "transect"), names(data))] <- 
    c("can", "tmnt", "tsct")
binwid <- 2.5
nbin <- 10 / binwid
binfac <- rep(1:nbin, each = 40 / nbin)
data$bin <- rep(binfac, length = nrow(data))
data$can <- as.character(ifelse(data$can == "barloy", 2, 1))
data$tmnt <- as.character(ifelse(data$tmnt == "control", 0, 1))
data$occ <- as.double(data$occ)
data <- aggregate(occ ~ can + bout + plot + tmnt + tsct + bin + cat, 
                  data, mean)
data$binid <- with(data, paste0(can, plot, tmnt, tsct, bin))
data$bout <- as.integer(data$bout)

# filter data ------------------------------------------------------------------
data <- # bapi data were collected with 10-cm rule
    with(data, data[!(can == 1 & cat == "bapi" & bout == 1), ])
# transects too close
data <- # pil p6 g t2 is too close
    with(data, data[!(can == 1 & plot == 6 & tmnt == 1 & tsct == 2), ])
data <- # bar p1 c t2 is too close
    with(data, data[!(can == 2 & plot == 1 & tmnt == 0 & tsct == 2), ])
data <- # bar p2 c t2 is too close
    with(data, data[!(can == 2 & plot == 2 & tmnt == 0 & tsct == 2), ])
data <- # bar p2 g t2 is too close
    with(data, data[!(can == 2 & plot == 2 & tmnt == 1 & tsct == 2), ])
data <- # bar p3 g t2 is too close
    with(data, data[!(can == 2 & plot == 3 & tmnt == 1 & tsct == 2), ])
data <- # bar p4 g t2 is too close
    with(data, data[!(can == 2 & plot == 4 & tmnt == 1 & tsct == 2), ])
# collected in different seasons
data <- # pilarcitos data were collected in winter
    with(data, data[!(can == 1 & cat != "bapi" & bout == 0), ])
data <- # barloy data were collected in winter
    with(data, data[!(can == 2 & cat != "bapi" & bout == 0), ])
data <- # barloy data do not have enough data from other bouts
    with(data, data[!(can == 2 & cat != "bapi" & bout == 1), ])
data <- # barloy data were collected in fall
    with(data, data[!(can == 2 & cat != "bapi" & bout == 4), ])

# export summary table ---------------------------------------------------------
mu <- aggregate(occ ~ cat + bout + tmnt, data = data, FUN = mean)
mu <- reshape(mu, # what data to reshape
              idvar = c("cat", "tmnt"), # remaining columns on left
              v.names = "occ", # names of new columns on right
              timevar = "bout", # values of new columns on right 
              direction = "wide") # going from long to wide
# give sortable column names
pre <- sort(unique(data$bout))
suf <- rep("_mu", length = length(pre))
colnames(mu)[-c(1:2)] <- paste0(pre, suf)

#define function for getting se for binomial distribution
get_se <- function(x) {
    k <- 1
    p <- mean(x)
    q <- 1 - p
    n <- length(x)
    se <- sqrt((p * q) / (k * n))
    return(se)
}
sig <- aggregate(occ ~ cat + bout + tmnt, data = data, FUN = get_se)
sig <- reshape(sig, 
               idvar = c("cat", "tmnt"), 
               v.names = "occ", 
               timevar = "bout", 
               direction = "wide")
pre <- sort(unique(data$bout))
suf <- rep("_sig", length = length(pre))
colnames(sig)[-c(1:2)] <- paste0(pre, suf)
summ_r <- round(cbind(mu[-(1:2)], sig[-(1:2)]) * 100)
summ_r <- summ_r[order(colnames(summ_r))]
summ_uns <- cbind(mu[1:2], summ_r)
# define row order
cats <- unique(data$cat)
rows_ctrl <- match(cats, summ_uns$cat)
rows_tmnt <- rows_ctrl + length(cats)
rows <- unlist(mapply(c, rows_tmnt, rows_ctrl, SIMPLIFY = F))
# order rows
summ <- summ_uns[rows, ]
# export table
write.csv(summ, paste0("./tables/summary.csv"), 
          row.names = F, quote = F)

# define functions
make_aic_table <- function(aic) {
    AIC <- aic$AIC
    aic_c <- AIC
    d_aic <- aic_c - min(aic_c)
    aic_w <- exp(-0.5 * d_aic) / sum(exp(-0.5 * d_aic))
    er <- max(aic_w) / aic_w
    ler <- log10(er)
    table <- cbind(aic, d_aic, aic_w, er, ler)
    table[order(table$AIC), ]
}
get_weighted_params <- function(aic_table, fits) {
    full <- fits[[match("mod_full", names(fits))]]
    fixef_names <- names(fixef(full))
    ranef_names <- c("binid_mu", "binid_sig")
    names <- c("cat", "model", fixef_names, ranef_names)
    params <- as.data.frame(matrix(NA, nrow = length(fits), ncol = length(names)))
    names(params) <- names
    fix_efs <- lapply(fits, fixef)
    ran_efs <- lapply(fits, function(x) ranef(x)[[1]][[1]])
    params$cat <- cat
    params$model <- names(fits)
    for (i in 1:length(fits)) {
        params[i, match(names(fix_efs[[i]]), names(params))] <- fix_efs[[i]]
        names(fix_efs[[1]])
        params[i, c("binid_mu", "binid_sig")] <- c(mean(ran_efs[[i]]), sd(ran_efs[[i]]))
    }
    return(params)
}
get_predicted <- function(aic_table, fits) {
    newdata <- data.frame("cat" = cat,
                          "binid" = id_fit,
                          "tmnt" = tmnt_fit,
                          "bout" = bout_fit,
                          "occ" = NA)
    occ_w <- 0
    for (i in 1:length(fits)) {
        occ <- predict(fits[[i]], newdata = newdata, type = "response",
                       allow.new.levels = TRUE) 
        occ <- occ * aic_table$aic_w[i]
        occ_w <- occ_w + occ 
    }
    newdata$occ <- occ_w
    return(newdata)
}

# define variables ------------------------------------------------
bouts <- sort(unique(data$bout))
smoothness <- 100
pad <- 1
bout_min <- min(data$bout) - pad
bout_max <- max(data$bout) + pad
bouts_seq <- seq(bout_min, bout_max, length = smoothness)
ids <- sort(unique(data$binid))
id_fit <- rep(ids, each = length(bouts_seq))
tmnt_fit <- as.character(ifelse(substr(id_fit, 3, 3) == 0, 0, 1))
bout_fit <- rep(bouts_seq, length(ids))

# run models 
order <- match(c("bapi", "bu", "ag", "onh", "onn", "inv"), cats)
cats <- cats[order] 
big_aic_table <- NULL
weighted_param_table <- NULL
newdata <- NULL
for (cat in cats) {
    data_sub <- data[data$cat == cat, ]
    mod_full  <- bglmer(occ ~ bout * tmnt + (1 | binid),
                    data = data_sub, fixef.prior = normal,
                    family = "binomial")
    mod_b_int  <- bglmer(occ ~ bout + bout:tmnt + (1 | binid),
                    data = data_sub, fixef.prior = normal,
                    family = "binomial")
    mod_t_b   <- bglmer(occ ~ tmnt + bout + (1 | binid),
                    data = data_sub, fixef.prior = normal,
                    family = "binomial")
    mod_t     <- bglmer(occ ~ tmnt + (1 | binid),
                    data = data_sub, fixef.prior = normal,
                    family = "binomial")
    mod_b     <- bglmer(occ ~ bout + (1 | binid),
                    data = data_sub, fixef.prior = normal,
                    family = "binomial")
    mod_null  <- bglmer(occ ~ 1 + (1 | binid),
                        data = data_sub, fixef.prior = NULL,
                        family = "binomial")
    aic <- AIC(mod_full, mod_b_int, mod_t_b, mod_t, mod_b, mod_null)
    aic_table <- make_aic_table(aic)
    aic_table$cat <- cat
    aic_table$model <- rownames(aic_table)
    aic_table <- aic_table[c("cat", "model", "df", "AIC", "d_aic", 
                            "aic_w", "er", "ler")]
    big_aic_table <- rbind(big_aic_table, aic_table)
    fits <- list(mod_full, mod_b_int, mod_t_b, mod_t, mod_b, mod_null)
    names(fits) <- c("mod_full", "mod_b_int", "mod_t_b", "mod_t", "mod_b", "mod_null")
    weighted_params <- get_weighted_params(aic_table, fits)
    weighted_param_table <- rbind(weighted_param_table, weighted_params)
    small_newdata <- 
        get_predicted(aic_table, fits)
    newdata <- rbind(newdata, small_newdata)
}
# write tables
write.csv(big_aic_table, "./tables/aic-table.csv", 
           row.names = F, quote = F)

# get weighted parameters
notcols <- match(c("binid_mu", "binid_sig"), names(weighted_param_table))
tab <- weighted_param_table[-notcols]
params <- reshape(tab, # data being reshaped
                idvar = c("cat", "model"), # repeats for each column being moved
                varying = 3:6, # the columns being moved
                times = names(tab)[3:6], # col names become values of resulting label col
                timevar = "param", # name of resulting label col
                v.names = "value", # name of resulting value col
                direction = "long")
rownames(params) <- NULL
luw <- big_aic_table[c(1, 2, 6)]
if(!"aic_w" %in% names(params))
    params <- merge(params, luw, by = c("cat", "model"))
params$valw <- params$value * params$aic_w
rowsri <- !is.na(params$value)
datri <- params[rowsri, ]
luri <- data.frame(with(datri, tapply(aic_w, list(cat, param), sum)))
luri$cat <- rownames(luri)
names(luri) <- c("cri", "bri", "xri", "tri", "cat")
luler <- luri 
colsler <- names(luri) != "cat"
luler[colsler] <- apply(luri[colsler], c(1, 2), function(x) log10(x / (1 - x)))
names(luler) <- c("cler", "bler", "xler", "tler", "cat")
rowsn <- params$cat == "bapi" & !is.na(params$value)
lun <- data.frame(table(params[rowsn, "param"]))
names(lun) <- c("param", "n")
if(!"n" %in% names(params))
    params <- merge(params, lun, by = "param")
params <- aggregate(valw ~ cat + param + n, params, sum)
params$val <- with(params, valw)
params <- params[-match(c("n", "valw"), names(params))]
params <- reshape(params, 
                  v.names = "val", 
                  idvar = "cat", 
                  timevar = "param",
                  direction = "wide")
names(params) <- gsub("val.", "", names(params))
names(params) <- c("cat", "x", "t", "b", "c") 
if(!"cri" %in% names(params)) {
    params <- merge(params, luri, by = "cat")
    params <- merge(params, luler, by = "cat")
}
params <- params[match(cats, params$cat), c("cat", sort(names(params)[-1])) ]
write.csv(params, "./tables/weighted-parameters.csv", 
           row.names = F, quote =F)


# plot -----------------------------------------------------------------
if (max(data$occ_perc <= 1)) {
    data$occ_perc <- data$occ * 100
    newdata$occ_perc <- newdata$occ * 100
}
cats <- c("bapi", "bu", "ag", NA, "onh", "onn")

png(paste0('./figures/model-fits.png'),
    width = 6, height = 6, units = "in", res = 300)
par(las = 1, mar = c(1, 1, 1, 0), oma = c(2, 3, 1, 1),
    mfcol = c(3, 2), xpd = NA, xaxs = "i", cex = 0.8)
xlab <- c("", "", "Grazing Bout", NA, "", "Grazing Bout")
ylab <- c("Percent Cover", "Percent Cover", "Percent Cover", NA, "", "")
xaxis <- c(FALSE, FALSE, TRUE, NA, FALSE, TRUE)
yaxis <- c(TRUE, TRUE, TRUE, NA, FALSE, FALSE) 
drawbox <- c(FALSE, FALSE, TRUE, NA, FALSE, FALSE)
wdat <- list(c(1:7), -c(1,5), -c(1, 5), NA, -c(1, 5), -c(1, 5))
labels <- c("Coyote Brush", "Bunchgrass", "Annual Grass",
            "Legend","Other Native Herbaceous", "Other Non-native")
drawleg <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
colt <- "#ff000088"
coltbox <- "#ff3333"
colc <- "#00000088"
colcbox <- "#333333"  
bwid <- 0.20
shift <- 0.12
lwd = 0.5
for (i in 1:length(cats)) {
    if (drawleg[i]) {
        plot(0, type = "n", axes = FALSE,
            xlim = c(bout_min, bout_max), ylim = c(0, 100), 
            xlab = "", ylab = "")
        lines(c(-0.7, 0.3) , c(60, 70), col = coltbox)
        lines(c(1, 2) , c(60, 70), col = colcbox)
        text(-.2, 85,  adj = c(0.5, 0.5),
            labels = "Treatment")
        text(1.5, 85,  adj = c(0.5, 0.5),
            labels = "Control")
        text(2.5, 70,  adj = c(0, 0.5),
            labels = "Regression Line")
        text(2.5, 58,  adj = c(0, 0.5),
            labels = expression(paste("for ", Bin[i])))
        lines(c(-0.2, -0.2) , c(10, 50), lwd = 0.5, col = coltbox)
        lines(c(1.5, 1.5) , c(10, 50), lwd = 0.5, col = colcbox)
        polygon(c(-0.3, -0.1, -0.1, -0.3), c(20, 20, 40, 40), 
                col = coltbox, border = coltbox)
        polygon(c(1.4, 1.6, 1.6, 1.4), c(20, 20, 40, 40), 
                col = colcbox, border = colcbox)
        text(2.5, 30, adj = c(0, 0.5),
             "Interquartile Range (box)\n& Total Range (whisker)\nof Groundcover")
        text(x = (bout_max - pad) / 2, 
             y = 105, 
             adj = c(0.5, 0),
             labels = labels[i])
    } else {
        cat <- cats[i]
        subc <- data[data$cat == cat & data$tmnt == 0, ]
        subt <- data[data$cat == cat & data$tmnt == 1, ]
        newc <- newdata[newdata$cat == cat & newdata$tmnt == 0, ]
        newt <- newdata[newdata$cat == cat & newdata$tmnt == 1, ]
        plot(0, type = "n", axes = F,
            xlim = c(bout_min, bout_max), ylim = c(0, 100), 
            xlab = "", ylab = "")
        title(xlab = xlab[i], line = 2)
        title(ylab = ylab[i], line = 3)
        boxplot(occ_perc ~ bout, 
                subt, 
                at = bouts[wdat[[i]]] - shift, 
                boxwex = bwid, 
                add = TRUE, 
                axes=FALSE, 
                col = coltbox,
                border = coltbox,
                staplewex = 0, 
                whisklty = "solid", 
                whisklwd = 0.5,
                range = 0,
                outline = F
        )
        boxplot(occ_perc ~ bout, 
                subc, 
                at = bouts[wdat[[i]]] + shift, 
                boxwex = bwid, 
                add = TRUE, 
                axes=FALSE, 
                col = colcbox,
                border = colcbox,
                staplewex = 0, 
                whisklty = "solid", 
                whisklwd = 0.5,
                range = 0,
                outline = F
        )
        for (bid in unique(newdata$binid)) {
            newcsub <- newc[newc$binid == bid, ]
            newtsub <- newt[newt$binid == bid, ]
            lines(newcsub$bout, newcsub$occ_perc, lwd = lwd, col = colc)
            lines(newtsub$bout, newtsub$occ_perc, lwd = lwd, col = colt)
        }
        if (yaxis[i]) {
            axis(2)
        }
        if (xaxis[i]) {
            axis(1, at = bout_min:bout_max, lwd.ticks = 0, labels = NA)
            axis(1, at = c(0, 1, 2, 3, 4, 5, 6))
        }
        if (drawbox[i]) box(bty = "l")
        text(x = (bout_max - pad) / 2, 
             y = 105, 
             adj = c(0.5, 0),
             labels = labels[i])
        if (cat == "bapi") {
            text(x = 0:6, y = rep(-5, 7), cex = 0.8,
                 c("P/B", "B", "P", "P", "B", "P", "P"))
        }
    }
}
dev.off()


