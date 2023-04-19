#!/bin/env Rscript

# Libraries
library(ggplot2)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
print(args)
tu_tsv = args[1]
wt_tsv = args[2]
out_png1 = args[3]
out_png2 = args[4]

tu_data <- read.table(tu_tsv, header=F)
wt_data <- read.table(wt_tsv, header=F)

data <- merge(tu_data, wt_data, by = "V2", suffixes = c("_1", "_2"))

# Plot for table 1
plot1 <- ggplot(data, aes(x = V2, y = V3_1)) +
  geom_point(color = "blue", alpha=0.1) +
  labs(title = "tumor", x = "X chromosome coordinates", y = "Read depth") +
  scale_y_continuous(limits = c(0, 100))

# Plot for table 2
plot2 <- ggplot(data, aes(x = V2, y = V3_2)) +
  geom_point(color = "red", alpha=0.1) +
  labs(title = "wild type", x = "X chromosome coordinates", y = "Read depth") +
  scale_y_continuous(limits = c(0, 100))

# Combine the plots using facet_wrap()
combined_plot <- plot1 + plot2 + facet_wrap(~ NULL, ncol=1, scales = "free_y")
# Display the plot
combined_plot
ggsave(out_png1, width=16, height=4, dpi=300)

data$log_fold <- log2(data$V3_1+1) - log2(data$V3_2+1)

ggplot(data, aes(x=V2, y=log_fold)) +
  geom_point() +
  ggtitle('Log2 read depth ratio')
ggsave(out_png2, width=16, height=4, dpi=300)

# ggplot(data, aes(x=V2, y=V3)) +
#   geom_point()  +
#   ggtitle(title)
# ggsave(out_png)


# ggplot(data, aes(x=V2, y=V3)) + 
#   geom_bar(stat = "identity") +
#   ggtitle(title) +
#   xlab('X chromosome coordinates') +
#   ylab('read depth')
# ggsave(out_png, width=12, height=4, dpi=300)






