library(tidyverse)
library(readxl)
library(ggthemes)
library(colorspace)
library(latex2exp)
library(scales)
library(directlabels)
library(geomtextpath)

library(showtext)
font_add_google("Noto Sans", "Noto Sans")

source("theme-mcm.R")

dat <- read_excel('./Cui2023TableS1.xlsx', sheet=1, range=cell_rows(3:175))

new_data <- data.frame(
  SiO2 = seq(45,65,by=1)
)

dat %>% ggplot(
  aes(x=SiO2, y=FeO)
) +
  geom_point() +
  geom_smooth(method='lm')


fit_FeO <- lm(FeO ~ SiO2, data=dat)
fit_Na2O <- lm(Na2O ~ SiO2, data=dat)
fit_CaO <- lm(CaO ~ SiO2, data=dat)
fit_MgO <- lm(MgO ~ SiO2, data=dat)
fit_Al2O3 <- lm(Al2O3 ~ SiO2, data=dat)
fit_TiO2 <- lm(TiO2 ~ SiO2, data=dat)
fit_K2O <- lm(K2O ~ SiO2, data=dat)
fit_P2O5 <- lm(P2O5 ~ SiO2, data=dat)
fit_MnO <- lm(MnO ~ SiO2, dat=dat)

new_data$FeO <- predict(fit_FeO, new_data)
new_data$Na2O <- predict(fit_Na2O, new_data)
new_data$CaO <- predict(fit_CaO, new_data)
new_data$MgO <- predict(fit_MgO, new_data)
new_data$Al2O3 <- predict(fit_Al2O3, new_data)
new_data$TiO2 <- predict(fit_TiO2, new_data)
new_data$K2O <- predict(fit_K2O, new_data)
new_data$P2O5 <- predict(fit_P2O5, new_data)
new_data$MnO <- predict(fit_MnO, new_data)

new_data <- new_data %>% mutate(Total = FeO+Na2O+CaO+SiO2+MgO+Al2O3+TiO2+K2O+P2O5+MnO)
new_data %>% write_delim("comps-interp.csv",delim=",")

