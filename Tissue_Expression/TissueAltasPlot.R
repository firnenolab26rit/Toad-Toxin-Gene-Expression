library(cowplot)
plot1 <- plot1 + theme(plot.margin = margin(15, 7, 7, 20))
plot2 <- plot2 + theme(plot.margin = margin(15, 7, 7, 20))
plot3 <- plot3 + theme(plot.margin = margin(30, 7, 7, 20))



# Stack plot1 and plot2 vertically
left_stack <- plot_grid(
  plot1, plot2,
  ncol = 1,
  align = "v",
  labels = c("A", "B"),     # Labels for the stacked plots
  label_size = 30,          # Font size
  label_fontface = "bold",   # Bold labels
  label_x = -0.005, 
  label_y = 1.0

)


# Combine left_stack with plot3 horizontally
final_plot <- plot_grid(
  left_stack, plot3,
  ncol = 2,
  rel_widths = c(2, 1),
  labels = c("", "C"),       # Label only for plot3, leave left_stack blank
  label_size = 30,
  label_fontface = "bold",
  label_y = c(0, 1.0)
)

# Display
final_plot + theme(plot.margin = margin(40, 10, 10, 10))

ggsave(filename = "Figure1TissueAltasFigure.jpeg", plot = final_plot, dpi = 600, width = 20, height = 8.27, units = "in")
ggsave(filename = "Figure1TissueAltasFigure.pdf", plot = final_plot, dpi = 600, width = 20, height = 8.27, units = "in", device = cairo_pdf)
ggsave(filename = "Figure1TissueAltasFigure.svg", plot = final_plot, dpi = 600, width = 20, height = 8.27, units = "in")

