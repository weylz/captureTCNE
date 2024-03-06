#!/usr/bin/env Rscript
# File    : fitTSD.R
# Time    : 2023/01/01 00:00:00
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : fit distribution of transcriptional signal


suppressMessages(library(fitdistrplus))
args <- commandArgs(TRUE)

message(
    cat("[", format(Sys.time(), "%a %b %d %I:%M:%S %p %Z %Y"), "] ", sep = ""),
    "Loading the distribution of transcription signal ...")
df <- read.table(args[1], sep = "\t", comment.char = "")
log_df <- log10(df[, 4])

fitn <- fitdist(log_df, "norm")
me <- round(as.numeric(fitn$estimate[1]), digits = 7)
sd <- round(as.numeric(fitn$estimate[2]), digits = 7)
# lower.tail: logical; if TRUE (default), probabilities are P[X <= x]
#           otherwise, P[X > x].
th <- round(
    qnorm(1 - as.numeric(args[2]), mean = me, sd = sd, lower.tail = FALSE),
    digits = 7)

write.table(round(10**th, digits = 7),
    "transcriptional.signal.cutoff.value.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)

pdf("transcriptional.signal.distribution.pdf", width = 5, height = 4)

hist(log_df, probability = TRUE, breaks = 100,
    col = "#3eb489", border = "#e5e4e2",
    xlab = "log10(TPM)", axes = TRUE,
    main = "Distribtion of Transcriptional Signal")
lines(density(log_df, bw = 0.2), col = "#9a9a9a", lwd = 3, lty = 1)
lines(density(rnorm(n = 1000000, mean = me, sd = sd), bw = 0.2),
    col = "#0038a8", lwd = 3, lty = 5)
lines(y = c(0, dnorm(th, mean = me, sd = sd)),
    x = c(th, th),
    col = "#ff1d52",
    lwd = 3)
legend("topright",
    c("empirical density curve", "fitted normal distribution"),
    col = c("#9a9a9a", "#0038a8"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")

invisible(dev.off())

message(
    cat("[", format(Sys.time(), "%a %b %d %I:%M:%S %p %Z %Y"), "] ", sep = ""),
    "Please check transcriptional.signal.distribution.pdf ",
    "for more information.")
