library(mgcv)
library(COVID19)
library(FactoMineR)
library(gplots)

setwd("./code")
OUTPUT_DIR <- "OUTPUT"

lag = 2

for (country in c("AUS", "BRA", "CHL", "COL", "CZE", "DEU", "JPN", "LTU", "ZAF")) {
    odir <- paste0("../output_interior-point_original_i1_d1_", country, sep = "")
    
    pref <- 1
    pname <- country
    
    sday <- as.Date("2021/4/1")
    eday <- as.Date("2021/11/10")
    
    n <- 3 # infected, dead, and recover

    A <- as.matrix(read.csv(paste0(odir, "/A1_", as.character(pref), ".csv"), header = FALSE))
    B <- as.matrix(read.csv(paste0(odir, "/B1_", as.character(pref), ".csv"), header = FALSE))
    K <- as.matrix(read.csv(paste0(odir, "/K1_", as.character(pref), ".csv"), header = FALSE))
    M <- as.matrix(read.csv(paste0(odir, "/MT_", as.character(pref), ".csv"), header = FALSE))
    MO <- as.matrix(read.csv(paste0(odir, "/MO_", as.character(pref), ".csv"), header = FALSE))

    nweek <- floor((as.numeric(eday - sday) + 1) / 7)
    nstep <- nweek - 1 # the estimation was conducted from the next step (week).

    X <- array(0, dim = c(n, nstep))
    X <- B %*% K %*% M[1:3, 1:nstep]
    X1 <- X[1, ]
    X2 <- X[2, ]
    X3 <- X[3, ]
    Y <- A %*% M[1:3, 1:nstep] + X
    YO <- MO[1:3, 2:(nstep + 1)]

    x <- covid19(country, level = 1, verbose = FALSE, amr = "../applemobilitytrends-2022-02-07.csv", start = "2021-03-31", end = "2021-11-22")
    xjpn <- covid19("JPN", level = 1, verbose = FALSE, amr = "../applemobilitytrends-2022-02-07.csv", start = "2021-03-31", end = "2021-11-22")
    
    nweek <- floor((as.numeric(eday - sday) + 1) / 7)
    nstep <- nweek - 1
    vnum <- 6 # covariates: week, people_fully_vaccinated, driving, walking, transit, hosp
    vnum <- vnum + 1

    xp1 <- x[, 27:30]
    xpjpn1 <- xjpn[, 27:30]

    xp0 <- data.frame(matrix(NA, ncol = 4, nrow = dim(x)[1]))
    names(xp0) <- names(x[, 27:30])
    for (i in 1:4) {
        if (max(xjpn[, 26 + i]) == min(xjpn[, 26 + i])) {
            xp0[, i] <- cut(x[, 26 + i], breaks = c(0, 50, 90, 100), labels = c("L", "M", "H"))
        } else {
            xp0[, i] <- cut(x[, 26 + i], quantile(c(0, xjpn[, 26 + i], 100), probs = seq(0, 1, 1 / 3)), labels = c("L", "M", "H"))
        }
    }

    xp3 <- MCA(xp0, ncp = 4, graph = FALSE)
    xp3_hcpc <- HCPC(xp3, nb.clust = 3, graph = FALSE)
    xp3_clust <- data.frame(date = x$date, xp3_hcpc$data.clust$clust)
    
    V <- array(0, dim = c(vnum, nstep))
    for (i in 1:nstep) {
        x1 <- x[as.Date(x$date) >= as.Date(sday + 7 * (i - lag)) & as.Date(x$date) < as.Date(sday + 7 * (i + 1 - lag)), ]
        V[2, i] <- sum(x1$people_fully_vaccinated, na.rm = T)
        V[3, i] <- sum(x1$driving, na.rm = T)
        V[4, i] <- sum(x1$walking, na.rm = T)
        V[5, i] <- sum(x1$transit, na.rm = T)
        V[6, i] <- sum(x1$hosp, na.rm = T)
        
        xp3w <- xp3_clust[as.Date(xp3_clust$date) >= as.Date(sday + 7 * (i - lag)) & as.Date(xp3_clust$date) < as.Date(sday + 7 * (i + 1 - lag)), ]
        xp3w1 <- table(xp3w[, 2])
        j <- which.max(xp3w1)
        V[7, i] <- names(xp3w1)[j]
        
    }
    V[1, ] <- 1:nstep
    test0 <- data.frame(
        week = as.numeric(V[1, ]), people_fully_vaccinated = as.numeric(V[2, ]),
        driving = as.numeric(V[3, ]), walking = as.numeric(V[4, ]), transit = as.numeric(V[5, ]), hosp = as.numeric(V[6, ]), policy = as.factor(V[7, ])
    )
    testall <- data.frame(X1 = X1, X2 = X2, X3 = X3, test0)
    vtext <- "-1 + week + people_fully_vaccinated + driving + walking + transit + hosp + policy "
    text <- paste("mod = gam(list(")
    for (i in 1:n) {
        text <- paste0(text, "X", as.character(i), " ~ ", vtext)
        if (i < n) {
            text <- paste0(text, ", ")
        }
    }
    text <- paste0(text, ") , family = mvn(d =", as.character(n), "), data = testall)")
    eval(parse(text = text))
    sm <- summary(mod)
    sink(file = paste0("../", OUTPUT_DIR, "/", country, "_lag", as.character(lag), ".txt"))
    print(sm)
    sink()
    sm$p.table[sm$p.table[, 4] >= 0.01, ] <- 0
    se <- sm$p.table[, "Estimate"]
    x <- rep(0, 27)
    names(x) <- c("week", "people_fully_vaccinated", "driving", "walking", "transit", "hosp", "policy1", "policy2", "policy3", "week.1", "people_fully_vaccinated.1", "driving.1", "walking.1", "transit.1", "hosp.1", "policy1.1", "policy2.1", "policy3.1", "week.2", "people_fully_vaccinated.2", "driving.2", "walking.2", "transit.2", "hosp.2", "policy1.2", "policy2.2", "policy3.2")
    x[names(se)] <- se
    y <- country
    assign(y, c(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27]))
}
z <- rbind(AUS, BRA, CHL, COL, CZE, DEU, JPN, LTU, ZAF)
saveRDS(z, paste0("../", OUTPUT_DIR, "/z_lag", as.character(lag), ".obj"))

## dendrogram (Figure 9)
z <- readRDS(paste0("../", OUTPUT_DIR, "/z_lag", as.character(lag), ".obj"))
resz <- PCA(z, ncp = 5, graph = FALSE)
test <- HCPC(resz, nb.clust = -1, graph = FALSE)
pdf(file = "./figure9.pdf")
plot(test, choice = "3D.map")
dev.off()

## heatmap (Figure 8)
z <- as.data.frame(readRDS(paste0("../", OUTPUT_DIR, "/z_lag", as.character(lag), ".obj")))[, 1:27]
z_sc <- as.matrix((scale(z)))
pdf(file = "./figure8.pdf")
thres  = qt((0.05/2), df=25)   # 31 - number of variables(6) - number of intercepts (0)
min_z_sc = min(z_sc)
max_z_sc = max(z_sc)
z_std = sd(as.matrix(z))
thres = thres * z_std
heatmap.2(z_sc,
    scale = "none", col = c("blue", "gray", "red"),
    breaks = c(min_z_sc, thres, (-thres), max_z_sc),
    trace = "none", density.info = "none",
    margin = c(10, 8),
    column_title = "Covariates", row_title = "Countries",
    srtCol=45
)
dev.off()
