# Create table of field data in long format
rm(list=ls())

# Import data
dir <- "./field-data/"
files <- file.path(dir, list.files(dir))
datl <- lapply(files, read.csv, stringsAsFactors = F, strip.white = T, 
               na.string = "")
dat <- do.call(rbind, datl)

# import species list
species <- read.csv('./species-list.csv', strip.white = T,
                    na.string = "", stringsAsFactors = F)
                   
codes <- species$code
cats <- species$cat

# check
obs <- unique(unlist(dat[6:17]))
obs[! obs %in% codes ]
dat[dat$shrub.1 == "ons/bapi" & !is.na(dat$shrub.1), ]

# Replace species with categories,
dat[6:17] <- sapply(dat[6:17], function(j) cats[match(j, codes)])

# Replace all occurrences with unique occurrences,
makeunique <- function (i, blank) {
    uni <- unique(as.character(i))
    blank[1:length(uni)] <- uni
    return(blank)
}
blank4 <- rep(NA, 4)
dat[, 6:9] <- t(apply(dat[, 6:9], 1, makeunique, blank4))
blank8 <- rep(NA, 8)
dat[, 10:17] <- t(apply(dat[, 10:17], 1, makeunique, blank8))

# one col per cat
cats <- sort(unique(cats))
wide <- dat[1:5]
wide[6:(5 + length(cats))] <- NA
names(wide)[6:ncol(wide)] <- cats

# Represent category occurrences with 1=occurrence and 0=nonoccurrence. 
codeoccs <- function (i) {
    cols <- match(i, cats)
    blank[cols] <- 1
    blank[is.na(blank)] <- 0
    return(blank)
}
blank <- rep(NA, length(cats))
wide[, 6:ncol(wide)] <- t(apply(dat[, 6:17], 1, codeoccs))

wide <- wide[c("canyon", "bout", "plot", "treatment", "transect", 
          "ag", "bapi", "bu", "inv", "onh", "onn")]
long <- reshape(wide, # what data to reshape
                direction = "long",
                varying = list(6:11),
                v.names = "occ",
                timevar = "cat",
                times = names(wide)[6:11])
long <- long[-match("id", names(long))]
write.csv(long, "data-long.csv", quote=F, row.names=F)
