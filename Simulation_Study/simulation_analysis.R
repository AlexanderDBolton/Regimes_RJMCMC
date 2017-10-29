setwd("~/Dropbox/Pairs_simulation_study")

actual_min_prop_similarities = numeric(1000)
mean_min_prop_similarities = numeric(1000)
quantiles_min_prop_similarities = list()
actual_similarities = numeric(1000)
mean_similarities = numeric(1000)
quantiles_similarities = list()
for (i in 1:1000) {
    filename = paste("comparison_results/comparison_results_", i - 1, ".txt", sep = "")
    X = read.table(filename)[,1]
    actual_min_prop_similarities[i] = X[5]
    mean_min_prop_similarities[i] = X[6]
    quantiles_min_prop_similarities[[i]] = X[7:11]
    actual_similarities[i] = X[12]
    mean_similarities[i] = X[13]
    quantiles_similarities[[i]] = X[14:18]
}

pdf("pairs_simulation_analysis.pdf")
par(mar = c(4.1, 4.1, 0.1, 0.1), cex = 2)
plot(actual_min_prop_similarities, mean_min_prop_similarities, xlab = "Actual similarity", ylab = "Mean estimated similarity", pch = ".")
abline(0, 1)
abline(0.0515, 1, lty = 2)
abline(-0.0515, 1, lty = 2)
dev.off()
