library('tidyverse')
library('ggthemes')

hil21 <- read_csv('./hillenbrand21.csv')

theme_mcm <- function() {
  font <- "Noto Sans"   #assign font family up front
  dgy <- "grey15"
  lgy <- "#cacaca"
  mrg <- 2
  theme_few() +
    theme(
      plot.title = element_text(
        family = font,
        size = 9,
        face = 'bold',
        margin = margin(b = mrg)
      ),
      
      plot.subtitle = element_text(family = font,
                                   size = 8),
      
      plot.caption = element_text(family = font,
                                  size = 8),
      
      panel.border = element_rect(
        color = dgy,
        size = 0.25,
        fill = NA
      ),
      
      
      axis.title = element_text(
        family = font,
        size = 9,
        color = dgy,
      ),
      axis.title.y.left = element_text(margin = margin(r = mrg)),
      axis.title.y.right = element_text(margin = margin(l = mrg)),
      axis.title.x.bottom =  element_text(margin = margin(t = mrg)),
      axis.title.x.top =  element_text(margin = margin(b = mrg)),
      
      axis.text = element_text(
        family = font,
        size = 8,
        color = dgy,
      ),
      
      strip.text = element_text(
        size = 8,
        face = 'bold',
        hjust = 0,
        margin = margin(b = mrg)
      ),
      
      legend.text = element_text(size = 8, color = dgy),
      
      legend.title = element_text(size = 9),
      
      panel.grid.major.y = element_line(
        color = lgy,
        size = 0.15,
        linetype = 1
      ),
      
      panel.grid.major.x = element_line(
        color = lgy,
        size = 0.15,
        linetype = 1
      ),
    )
}

# Set theme globally
theme_set(theme_mcm())

car22 <- data.frame(
  age=c(0, 
        11.4,
        11,
        11.5,
        39.3),
  age_sd=c(0.7,
           0.5,
           0.6,
           3.3),
  thickness=c(55, 
              38.3,
              38.9,
              42.7,
              34.8),
  thickness_min=c(0,
                  33.5,
                  37.4,
                  34,
                  30.4),
  thickness_max=c(0,
                  43.8,
                  43.1,
                  48,
                  39.2),
  name='Acongagua'
) %>% as.tibble

tMyr <- seq(0,50)
hkm <- seq(80,30)

car22 %>%
  ggplot(aes(x=age,y=thickness,colour=name)) + 
    geom_point() +
    geom_abline(intercept=80,slope=-1) +
    labs(x="Age (Ma)", y="Thickness (km)") +
    coord_cartesian(xlim = c(0, 50), ylim = c(30, 80)) +
    theme_mcm()

hil21 %>%
  ggplot(aes(x=age,y=thickness,colour=name)) + 
  geom_point() +
  geom_abline(intercept=80,slope=-1) +
  labs(x="Age (Ma)", y="Thickness (km)") +
  coord_cartesian(xlim = c(0, 50), ylim = c(30, 80)) +
  theme_mcm()
  
  