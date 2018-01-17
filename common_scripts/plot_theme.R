theme_pub = function(base_size = 15, font = 'Helvetica') {
  txt = element_text(size = base_size, colour = 'black', face = 'plain')
  bold_txt = element_text(size = base_size, colour = 'black', face = 'bold')
  
  theme_classic(base_size = base_size, base_family = font)  +
    theme(
      legend.key = element_blank(),
      strip.background = element_blank(),
      text = txt,
      plot.title = txt,
      axis.title = bold_txt,
      axis.text = txt,
      legend.title = bold_txt,
      legend.text = txt,
      legend.position = 'bottom')
}
  