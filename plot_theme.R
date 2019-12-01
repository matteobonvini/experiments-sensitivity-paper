require(ggplot2)

if (.Platform$OS.type != "unix"){
  windowsFonts(Times = windowsFont("TT Times New Roman"))
}

our_theme <- theme(plot.title = element_text(family = "Times", 
                                             size = 18,
                                             color = "black",
                                             hjust = 0.5),
                   panel.background = element_rect(fill = "white"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.x = element_text(family = "Times", 
                                              size = 12,
                                              color = "black",
                                              hjust = 0.5),
                   axis.text.y = element_text(family = "Times", 
                                              size = 12,
                                              color = "black"),
                   axis.title = element_text(family = "Times", 
                                             size = 14),
                   legend.title = element_text(family = "Times", 
                                               size = 14),
                   legend.title.align = 0.5,
                   legend.text = element_text(family = "Times", 
                                              size = 12,
                                              color = "black"),
                   legend.background = element_blank(),
                   legend.key = element_rect(fill = NA, colour = NA, size = 0.25),
                   plot.caption = element_text(family = "Times", 
                                               size = 14))
