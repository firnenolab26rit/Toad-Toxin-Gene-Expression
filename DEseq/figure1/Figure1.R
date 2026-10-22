# Figure 1 Plot

plot1 <- read_rds("../WithinTissue/WithinTissue_LinePlot.rds")

plot2 <- read_rds("../BetweenTissue/BetweenTissue_BarPlot.rds")
plot2 <- plot2 + labs(x = "Time Points")


plot1 + theme(plot.margin = margin(40, 10, 10, 10))
plot2 + theme(plot.margin = margin(40, 10, 10, 10))

Figure1 <- plot_grid(
  plot1,
  plot2,
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 30,
  hjust = -0.2,
  vjust = 1.1
)

ggsave(filename = "Figure2DE.jpeg", plot = Figure1, dpi = 600, width = 20, height = 8.27, units = "in")
ggsave(filename = "Figure2DE.pdf", plot = Figure1, dpi = 600, width = 20, height = 8.27, units = "in", device = cairo_pdf)
ggsave(filename = "Figure2DE.svg", plot = Figure1, dpi = 600, width = 20, height = 8.27, units = "in")



#write_csv(plot1$data, file = "WithinTissueDEGs_WaldResults.csv")
#write_csv(plot2$data, file = "BetweenTissueDEGs_WaldResults.csv")
