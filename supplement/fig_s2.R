#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')

#Load data
linear_fit = read.csv( 'raw_data/fitted_ercc.csv', header = TRUE, strip.white=TRUE)
measured_vals = read.csv('raw_data/obs_ercc.csv', header = TRUE, strip.white=TRUE)

#Plot data
plot = ggplot() +
  geom_point(data = measured_vals, aes(x=log2(x), y=y), size = 1.5, alpha = 0.7) +
  geom_line(data = linear_fit, aes(x=log2(x), y=y), size = 1.2) +
  theme_pub() +
  labs(x = 'log2(estimated count)', y = 'C_t value')

pdf('supplement/s5/fig_s5.pdf', useDingbats = FALSE, width = 2, height = 2)
print(plot)
dev.off()
