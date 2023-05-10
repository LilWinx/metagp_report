library("tidyverse")
library("ggplot2")
require("scales")

cov <- read.delim("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/images/21-P0542-M015.depth.txt")
data = subset(cov, select = -c(ID))
plot <- ggplot(data, aes(x=Position, y=Coverage)) +
  geom_area(fill="skyblue") +
  labs(
    title = "Test",
    x = "Base position in genome",
    y = "Fold10 Coverage"
  ) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 8)
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(trans='log10') +
  theme_bw()

plot
ggsave(
  plot = plot, 
  file = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/images/21-P0542-M015.png", 
  width = 1500, 
  height = 400,
  units = "px",
  dpi = 300,
  scale = 2
  )