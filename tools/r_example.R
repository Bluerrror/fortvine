library(fortvine)
x <- diff(log(EuStockMarkets))
u <- pseudo_obs(x)
pair <- bicop(u[, c("DAX", "CAC")], family_set = "all", cores = 2)
fit <- vinecop(u, family_set = "all", cores = 2)
set.seed(1)
sim <- rvinecop(nrow(u), fit, cores = 2)
w <- rosenblatt(u, fit, cores = 2)

blue <- "#1F5FA6"; green <- "#2E8B57"; gold <- "#C7922B"
ramp <- colorRampPalette(c("#FFFFFF", "#E3EDF8", "#AFCBEA", "#5F95D0", "#1F5FA6", "#0D2C52"))(64)
z <- seq(-3.2, 3.2, length.out = 161)
zz <- as.matrix(expand.grid(z, z))
dens <- matrix(dbicop(pnorm(zz), pair) * dnorm(zz[, 1]) * dnorm(zz[, 2]), length(z))
q <- function(v) qnorm(v)

png("docs/assets/r_example.png", width = 7, height = 2.6, units = "in", res = 140,
    type = "cairo", family = "Times New Roman", bg = "white")
par(mfrow = c(1, 3), mar = c(3.4, 3.4, 2.2, 1.2), mgp = c(2.1, 0.6, 0), cex = 0.78,
    las = 1, tcl = -0.3)
image(z, z, dens, col = ramp, xlab = "DAX (normal scores)", ylab = "CAC (normal scores)",
      xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2), useRaster = TRUE)
contour(z, z, dens, add = TRUE, drawlabels = FALSE, col = "#0D2C52", lwd = 0.5,
        levels = c(0.005, 0.02, 0.05, 0.1, 0.2, 0.3))
box()
title(main = sprintf("(a) selected pair copula: %s %d", toupper(pair$family), pair$rotation),
      font.main = 1, cex.main = 1)
plot(q(u[, "DAX"]), q(u[, "FTSE"]), pch = 16, cex = 0.3, col = adjustcolor(gold, 0.55),
     xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), xlab = "DAX (normal scores)",
     ylab = "FTSE (normal scores)")
points(q(sim[, 1]), q(sim[, 4]), pch = 16, cex = 0.3, col = adjustcolor(blue, 0.45))
legend("topleft", c("observed", "simulated from the vine"), pch = 16, pt.cex = 0.8,
       col = c(gold, blue), bty = "n", cex = 0.9)
title(main = "(b) data and vine simulation", font.main = 1, cex.main = 1)
plot(w[, 1], w[, 4], pch = 16, cex = 0.3, col = adjustcolor(green, 0.55), xlim = c(0, 1),
     ylim = c(0, 1), xlab = expression(w[1]), ylab = expression(w[4]), xaxs = "i", yaxs = "i")
tau <- cor(w[, 1], w[, 4], method = "kendall")
title(main = bquote("(c) Rosenblatt transform, " * tau == .(sprintf("%.3f", tau))), font.main = 1, cex.main = 1)
dev.off()
cat(cor(w[, 1], w[, 4], method = "kendall"), "\n")
