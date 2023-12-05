library(tidyverse)
library(ggthemes)
library(colorspace)
library(latex2exp)
library(scales)
library(directlabels)
library(geomtextpath)

library(showtext)
font_add_google("Noto Sans", "Noto Sans")

source('theme-mcm.R')


T0string = "$T_{0}$ (°C)"

bous_plag_out <- data.frame(
  plag_out_pressure = c(1.1036585365853657, 2.2012195121951224),
  plag_out_temperature = c(521.0053859964092, 1000.3590664272889)
) %>%
  as_tibble

hacker_2003 <- data.frame(
  plag_out_pressure = c(#  4.952095808383233,
    #  1.724550898203593,
    1.1158536585365857,
    1.451219512195122,
    1.8353658536585367),
  plag_out_temperature = c(#  403.1746031746031,
    #  489.94708994708986,
    521.0053859964092,
    670.7360861759424,
    1001.4362657091561)
) %>%
  as_tibble

bous_eclogite <- data.frame(
  plag_out_pressure = c(#   2.562874251497006,
    #  1.7664670658682633,
    1.25,
    2.353658536585366),
  plag_out_temperature = c(#   478.3068783068782,
    #   480.42328042328035,
    515.6193895870736
    , 1003.5906642728902)
) %>%
  as_tibble

dat <- read_csv('./_critical.csv') %>%
  mutate(critical_depth = ifelse(critical_depth > 80e3, NA, critical_depth)) %>%
  mutate(S0 = 6000) %>%
  mutate(h0 = 50e3) %>%
  mutate(v0 = 1e-3 / 3.154e7) %>%
  mutate(t0 = h0 / v0) %>%
  mutate(r0 = Da / S0 / h0 * v0 * rho0)


dat %>%
  group_by(composition) %>%
  mutate(mean_depth = mean(critical_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(composition) %>%
  ggplot(aes(
    x = Da,
    y = critical_depth / 1e3,
    colour = T0 - 273.15,
    group = interaction(T0, composition)
  )) +
  scale_colour_continuous_diverging(
    palette = "Blue-Red 3",
    l1 = 50,
    l2 = 90,
    p1 = 0.5,
    p2 = 1.2,
    mid = 550
  ) +
  geom_line(linewidth = 0.3, alpha = 1) +
  geom_point(size = 0.4) +
  scale_y_reverse(limits = c(85, 25), breaks = seq(30, 80, by = 10)) +
  scale_x_continuous(
    trans = "log10",
    breaks = trans_breaks('log10', function(x)
      10 ^ x),
    labels = trans_format('log10', math_format(10 ^ .x))
  ) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  labs(
    y = TeX("stable crustal thickness (km)"),
    x = TeX("Damköhler number"),
    colour = TeX(T0string)
  ) +
  facet_wrap( ~ composition)

ggsave('plots/critical_by_composition.pdf', device = cairo_pdf)
ggsave('plots/critical_by_composition.png')

dat %>%
  group_by(composition) %>%
  mutate(mean_depth = mean(critical_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(composition) %>%
  ggplot(aes(
    x = T0 - 273.15,
    y = critical_depth / 1e3,
    colour = log10(Da),
    group = interaction(Da, composition)
  )) +
  scale_colour_viridis_c() +
  geom_line(linewidth = 0.3, alpha = 1) +
  geom_point(size = 0.4) +
  scale_y_reverse(limits = c(80, 30), breaks = seq(30, 80, by = 10)) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  labs(
    y = TeX("stable crustal thickness (km)"),
    x = TeX(T0string),
    colour = TeX("log Da")
  ) +
  facet_wrap( ~ composition)

ggsave('plots/critical_by_composition_alt.pdf', device = cairo_pdf)
ggsave('plots/critical_by_composition_alt.png')

mafic <- 'mackwell_1998_maryland_diabase_norm'
felsic <- 'sammon_2021_lower_crust_norm'

c1low <-
  dat %>% filter(composition == mafic) %>% filter(T0 == min(T0)) %>% arrange(Da)
c1hi <-
  dat %>% filter(composition == mafic) %>% filter(T0 == max(T0)) %>% arrange(-Da)
c1x <- append(c1low$Da, c1hi$Da)
c1y <- append(c1low$critical_depth, c1hi$critical_depth)

c2low <-
  dat %>% filter(composition == felsic) %>% filter(T0 == min(T0)) %>% arrange(Da)
c2hi <-
  dat %>% filter(composition == felsic) %>% filter(T0 == max(T0)) %>% arrange(-Da)
c2x <- append(c2low$Da, c2hi$Da)
c2y <- append(c2low$critical_depth, c2hi$critical_depth)

c1y[is.na(c1y)] <-  85.e3
c2y[is.na(c2y)] <- 85.e3

ggplot() +
  geom_polygon(aes(x = c2x, y = c2y / 1e3, fill = 'Lower crust'), alpha =
                 1) +
  geom_polygon(aes(x = c1x, y = c1y / 1e3, fill = 'Maryland Diabase'), alpha =
                 1) +
  scale_fill_manual(values = c("#ffeeee", "#dfeeff")) +
  geom_textline(
    data = dat %>% filter(composition == c1low$composition),
    aes(
      x = Da,
      y = critical_depth / 1e3,
      group = T0,
      colour = '2cyan',
      label = round(T0 - 273.15),
      alpha = T0
    ),
    hjust = 0.75,
    size = 3,
    show.legend = FALSE,
    text_smoothing = 30
  ) +
  geom_textline(
    data = dat %>% filter(composition == c2low$composition),
    aes(
      x = Da,
      y = critical_depth / 1e3,
      group = T0,
      colour = '1red',
      label = round(T0 - 273.15),
      alpha = T0
    ),
    hjust = 0.9,
    size = 3,
    show.legend = FALSE,
    text_smoothing = 30
  ) +
  scale_colour_manual(values = c("#ff4444", "#00aacc")) +
  scale_alpha(range = c(0.4, 1.0)) +
  scale_y_reverse(breaks = seq(30, 80, by = 10), expand = c(0, 0)) +
  scale_x_continuous(
    trans = "log10",
    breaks = trans_breaks('log10', function(x)
      10 ^ x),
    labels = trans_format('log10', math_format(10 ^ .x)),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(1e-3, 1e6), ylim = c(80, 30)) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  labs(
    y = TeX("Stable depth (km)"),
    x = TeX("Damköhler number"),
    colour = TeX(T0string),
    fill = 'Composition'
  ) +
  theme(legend.position = "bottom")

ggsave(
  'plots/critical_shaded.pdf',
  device = cairo_pdf,
  width = 6,
  height = 4
)
ggsave('plots/critical_shaded.png',
       width = 6,
       height = 4)

c1low <- dat %>%
  filter(composition == mafic) %>%
  filter(Da == min(Da)) %>%
  arrange(T0)

c1hi <- dat %>%
  filter(composition == mafic) %>%
  filter(Da == max(Da)) %>%
  arrange(-T0)

c1x <- append(c1low$T0 - 273.15, c1hi$T0 - 273.15)
c1y <- append(c1low$critical_depth, c1hi$critical_depth)

c2low <- dat %>%
  filter(composition == felsic) %>%
  filter(Da == min(Da)) %>%
  arrange(T0)
c2hi <- dat %>%
  filter(composition == felsic) %>%
  filter(Da == max(Da)) %>%
  arrange(-T0)

c2x <- append(c2low$T0 - 273.15, c2hi$T0 - 273.15)
c2y <- append(c2low$critical_depth, c2hi$critical_depth)

c1y[is.na(c1y)] <-  95.e3
c2y[is.na(c2y)] <- 95.e3

ggplot() +
  geom_polygon(aes(x = c1x, y = c1y / 1e3, fill = "Maryland diabase"), alpha =
                 1) +
  geom_polygon(aes(x = c2x, y = c2y / 1e3, fill = "Int.lower crust"), alpha =
                 1) +
  scale_fill_manual(values = c("#ffeeee", "#dfeeff")) +
  geom_textline(
    data = dat %>% filter(composition == mafic &
                            Da <= 300 &
                            (Da * 1e9) %% 3 != 0) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), NaN, critical_depth
                            )),
    aes(
      x = T0 - 273.15,
      y = critical_depth / 1e3,
      group = Da,
      colour = "Maryland diabase",
      label = Da,
      #alpha = log10(Da)
    ),
    hjust = 0.3,
    size = 2.5,
    show.legend = FALSE,
    text_smoothing = 30,
    linewidth = 0.4
  ) +
  geom_line(
    data = dat %>% filter(composition == mafic &
                            Da <= 300 &
                            (Da * 1e9) %% 3 == 0) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), NaN, critical_depth
                            )),
    aes(
      x = T0 - 273.15,
      y = critical_depth / 1e3,
      group = Da,
      colour = "Maryland diabase",
      #alpha = log10(Da)
    ),
    linewidth = 0.4,
    linetype = "32",
    show.legend = FALSE,
    alpha = 0.5
  ) +
  geom_textline(
    data = dat %>% filter(composition == felsic &
                            Da <= 300 &
                            (Da * 1e9) %% 3 != 0) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), NaN, critical_depth
                            )),
    aes(
      x = T0 - 273.15,
      y = critical_depth / 1e3,
      group = Da,
      colour = "Int.lower crust",
      label = Da,
      alpha = log10(Da)
    ),
    hjust = 0.45,
    size = 2.5,
    show.legend = FALSE,
    text_smoothing = 30,
    linewidth = 0.4
  ) +
  geom_line(
    data = dat %>% filter(composition == felsic &
                            Da <= 300 &
                            (Da * 1e9) %% 3 == 0) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), NaN, critical_depth
                            )),
    aes(
      x = T0 - 273.15,
      y = critical_depth / 1e3,
      group = Da,
      colour = "Int.lower crust",
      label = Da,
      #alpha = log10(Da)
    ),
    show.legend = FALSE,
    linewidth = 0.4,
    linetype = "32",
    alpha = 0.5
  ) +
  geom_line(
    data = dat %>% filter(composition == mafic &
                            Da == 1e3) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), 85.e3, critical_depth
                            )),
    aes(x = T0 - 273.15, y = critical_depth / 1e3),
    colour = "#004466"
  ) +
  geom_line(
    data = dat %>% filter(composition == felsic &
                            Da == 1e3) %>% mutate(critical_depth = ifelse(
                              is.na(critical_depth), 85.e3, critical_depth
                            )),
    aes(x = T0 - 273.15, y = critical_depth / 1e3),
    colour = "#aa0000"
  ) +
  scale_colour_manual(values = c("#ff4444", "#00aacc")) +
  scale_alpha(range = c(0.3, 1.0)) +
  scale_y_reverse(breaks = seq(30, 80, by = 10), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 50),
    expand = c(0, 0),
    labels = ifelse(
      seq(0, 1000, by = 50) %% 100 == 0,
      format(seq(0, 1000, by = 50), digits = 1, nsmall = 0),
      ""
    )
  ) +
  coord_cartesian(ylim = c(80, 30)) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.05, "cm"),
    mid = unit(0.05, "cm"),
    long = unit(0.1, "cm")
  ) +
  labs(
    y = TeX("Stable depth (km)"),
    x = TeX(T0string),
    colour = TeX("log Da"),
    fill = 'Composition'
  ) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank()) +
  theme(legend.position = c(0.15, 0.15))

ggsave(
  'plots/critical_shaded_alt.pdf',
  device = cairo_pdf,
  width = 3.5,
  height = 3.1
)
ggsave('plots/critical_shaded_alt.png',
       width = 3.5,
       height = 3.1)

dat %>%
  mutate(critical_depth = ifelse(critical_depth > 80e3, NA, critical_depth)) %>%
  group_by(composition) %>%
  mutate(mean_depth = mean(critical_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(composition) %>%
  ggplot(aes(
    x = Da,
    y = plag_out_depth / 1e3,
    colour = T0 - 273.15,
    group = interaction(T0, composition)
  )) +
  scale_colour_continuous_diverging(
    palette = "Blue-Red 3",
    l1 = 50,
    l2 = 90,
    p1 = 0.5,
    p2 = 1.2,
    mid = 550
  ) +
  geom_line(linewidth = 0.3, alpha = 1) +
  geom_point(size = 0.4) +
  scale_y_reverse(limits = c(85, 25), breaks = seq(30, 80, by = 10)) +
  scale_x_continuous(
    trans = "log10",
    breaks = trans_breaks('log10', function(x)
      10 ^ x),
    labels = trans_format('log10', math_format(10 ^ .x))
  ) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  labs(
    y = TeX("Plag-out depth (km)"),
    x = TeX('Damköhler number'),
    colour = TeX(T0string)
  ) +
  facet_wrap( ~ composition)

ggsave('plots/plag_out_depth_by_composition.pdf', device = cairo_pdf)
ggsave('plots/plag_out_depth_by_composition.png')

dat %>%
  mutate(critical_depth = ifelse(critical_depth > 80e3, NA, critical_depth)) %>%
  group_by(composition) %>%
  mutate(mean_depth = mean(critical_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(composition) %>%
  ggplot(aes(
    x = T0 - 273.15,
    y = plag_out_depth / 1e3,
    colour = log10(Da),
    group = interaction(Da, composition)
  )) +
  scale_color_viridis_c() +
  geom_line(linewidth = 0.3, alpha = 1) +
  geom_point(size = 0.4) +
  scale_y_reverse(limits = c(85, 25), breaks = seq(30, 80, by = 10)) +
  annotation_logticks(
    sides = "b",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  labs(
    y = TeX("Plag-out depth (km)"),
    x = TeX(T0string),
    colour = TeX("log Da")
  ) +
  facet_wrap( ~ composition)

ggsave('plots/plag_out_depth_by_composition_alt.pdf', device = cairo_pdf)
ggsave('plots/plag_out_depth_by_composition_alt.png')

base_plot <- dat %>%
  filter(Da < 3000) %>%
  ggplot(
    aes(
      x = plag_out_temperature - 273.15,
      y = plag_out_pressure / 1e4,
      colour = as.factor(Da)
    ),
    group = "1"
  ) +
  geom_line(
    data = hacker_2003,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Hacker 2003 Basalt"),
    colour = "darkred",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_line(
    data = bous_plag_out,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Bousquet plag-out"),
    colour = "red",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_line(
    data = bous_eclogite,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Bousquet eclogite"),
    colour = "orange",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_textline(
    aes(label = Da),
    size = 3,
    show.legend = FALSE,
    text_smoothing = 30
  ) +
  #geom_dl(aes(label=factor(Da)), method=list("last.points",cex=0.6)) +
  scale_colour_grey(start = 0.7, end = 0, name = "Da") +
  labs(y = TeX("Plag-out pressure (GPa)"),
       x = TeX("Temperature (°C)")) +
  scale_y_continuous(breaks = seq(0, 3, by = .5), expand = c(0, 0), ) +
  scale_x_continuous(breaks = seq(500, 900, by = 100), expand = c(0, 0)) +
  coord_cartesian(xlim = c(400, 1000), ylim = c(1, 2.5)) +
  theme(aspect.ratio = 1) +
  scale_size_manual(" ",
                    values = rep(1, 3),
                    guide = guide_legend(override.aes = list(colour = c(
                      "orange", "red", "darkred"
                    ))))

base_plot + facet_wrap( ~ composition)

ggsave(
  'plots/plag_out_pt.pdf',
  device = cairo_pdf,
  width = 12,
  height = 12
)
ggsave('plots/plag_out_pt.png', width = 12, height = 12)

base_plot %+%
  filter(
    dat,
    Da < 1000 &
      (
        composition == 'mackwell_1998_maryland_diabase_norm' |
          composition == 'sammon_2021_lower_crust_norm'
      )
  ) +
  facet_wrap( ~ composition) +
  theme(legend.position = "bottom")

ggsave(
  'plots/plag_out_pt_MD_LC.pdf',
  device = cairo_pdf,
  width = 5,
  height = 4
)
ggsave('plots/plag_out_pt_MD_LC.png',
       width = 5,
       height = 4)

base_plot <- dat %>%
  mutate(critical_pressure = ifelse(critical_pressure > 28e3, NaN, critical_pressure)) %>%
  ggplot(
    aes(
      x = critical_temperature - 273.15,
      y = critical_pressure / 1e4,
      group = interaction(Da, composition),
      colour = as.factor(Da)
    )
  ) +
  geom_line(
    data = hacker_2003,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Hacker 2003 Basalt"),
    colour = "darkred",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_line(
    data = bous_plag_out,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Bousquet plag-out"),
    colour = "red",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_line(
    data = bous_eclogite,
    aes(x = plag_out_temperature, y = plag_out_pressure, size = "Bousquet eclogite"),
    colour = "orange",
    linetype = '42',
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_textline(
    aes(label = Da),
    size = 3,
    show.legend = FALSE,
    text_smoothing = 30
  ) +
  geom_textline(aes(label = Da), size = 3, show.legend = FALSE) +
  scale_colour_grey(start = 0.8, end = 0, name = "Da") +
  scale_y_continuous(breaks = seq(0, 3, by = .5), expand = c(0, 0), ) +
  scale_x_continuous(breaks = seq(500, 900, by = 200), expand = c(0, 0)) +
  coord_cartesian(xlim = c(400, 1000), ylim = c(0.5, 2.5)) +
  labs(y = TeX("Critical pressure (GPa)"),
       x = TeX("$Temperature$ (°C)")) +
  theme(aspect.ratio = 1) +
  scale_size_manual(" ",
                    values = rep(1, 3),
                    guide = guide_legend(override.aes = list(colour = c(
                      "orange", "red", "darkred"
                    ))))

base_plot + facet_wrap( ~ composition)

ggsave(
  'plots/critical_pt.pdf',
  device = cairo_pdf,
  width = 12,
  height = 12
)
ggsave('plots/critical_pt.png', width = 12, height = 12)


base_plot %+%
  filter(
    dat,
    Da %in% c(0.001, 0.01, 0.1, 0.3, 1, 3, 10, 30, 100, 1000) &
      (
        composition %in% c(
          'mackwell_1998_maryland_diabase_norm',
          'sammon_2021_lower_crust_norm',
          'hacker_2015_md_xenolith_norm'
        )
      ) & (critical_pressure < 30e3)
  ) +
  facet_wrap(
    ~ factor(
      composition,
      levels = c(
        'mackwell_1998_maryland_diabase_norm',
        'hacker_2015_md_xenolith_norm',
        'sammon_2021_lower_crust_norm'
      )
    ),
    labeller = function(labels)
      list(LETTERS[1:nrow(labels)])
  ) +
  theme(legend.position = "bottom")

ggsave('plots/critical_pt_MD_LC.pdf',
       device = cairo_pdf,
       width = 6.5)
ggsave('plots/critical_pt_MD_LC.png', width = 6.5)

dat %>%
  filter(Da %in% c(0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6)) %>%
  filter(
    composition %in% c(
      "hacker_2015_md_xenolith_norm",
      "mackwell_1998_maryland_diabase_norm",
      "sammon_2021_lower_crust_norm",
      "sammon_2021_deep_crust_norm"
    )
  ) %>%
  mutate(name = composition) %>%
  mutate(name = str_replace(name, "sammon_2021_deep_crust_norm", "57.6")) %>%
  mutate(name = str_replace(name, "sammon_2021_lower_crust_norm", "53.3")) %>%
  mutate(name = str_replace(name, "hacker_2015_md_xenolith_norm", "52.1")) %>%
  mutate(name = str_replace(name, "mackwell_1998_maryland_diabase_norm", "51.6")) %>%
  ggplot(
    aes(
      x = (Da),
      y = max_densification_rate_10kyr,
      ymax = max_densification_rate_10kyr,
      ymin = max_densification_rate_10kyr,
      fill = factor(name),
      colour = (T0 - 273.15)
    )
  ) +
  annotate(
    "text",
    x = .03,
    y = 2.e2,
    label = "fluid-poor",
    colour = "#6699ff"
  ) +
  annotate(
    "segment",
    x = 0.003,
    xend = 0.25,
    y = 1.8e2,
    yend = 1.8e2,
    colour = '#aabbff',
    size = 0.5,
    arrow = arrow(
      type = "closed",
      length = unit(0.01, "npc"),
      ends = "both"
    )
  ) +
  annotate(
    "text",
    x = 1,
    y = 2.e2,
    label = "~1% H2O",
    colour = "#ff6666"
  ) +
  annotate(
    "segment",
    x = 0.35,
    xend = 3,
    y = 1.8e2,
    yend = 1.8e2,
    colour = '#ffaaaa',
    size = 0.5,
    arrow = arrow(
      type = "closed",
      length = unit(0.01, "npc"),
      ends = "both"
    )
  ) +
  annotate(
    "text",
    x = 3e5,
    y = 2e2,
    label = "fluid-rich",
    colour = "#66aa88"
  ) +
  annotate(
    "segment",
    x = 3e4,
    xend = 3e6,
    y = 1.8e2,
    yend = 1.8e2,
    colour = '#77aa99',
    size = 0.5,
    arrow = arrow(
      type = "closed",
      length = unit(0.01, "npc"),
      ends = "both"
    )
  ) +
  geom_abline(
    intercept = log10(40),
    slope = 0,
    size = 0.6,
    alpha = 0.5,
    linetype = "longdash",
    colour = 'black'
  ) +
  geom_abline(
    intercept = log10(80),
    slope = 0,
    size = 0.6,
    alpha = 0.5,
    linetype = "longdash",
    colour = 'black'
  ) +
  scale_y_continuous(
    limits = c(1, 2.e2),
    trans = "log10",
    breaks = c(1e0, 1e1, 1e2, 1e3, 1e4),
    labels = trans_format('log10', math_format(10 ^ .x))
  ) +
  scale_x_log10(
    expand = expansion(add = .5),
    breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
    labels = trans_format('log10', math_format(10 ^ .x))
  ) +
  
  annotation_logticks(
    sides = "l",
    size = 0.25,
    short = unit(0.1, "cm"),
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm")
  ) +
  stat_summary(
    fun.max = max,
    fun.min = min,
    fun = mean,
    position = position_dodge(0.75),
    linewidth = 4,
    size = 0,
    colour = "#eaeaea",
    lineend = 'butt',
    show.legend = FALSE
  ) +
  geom_point(position = position_jitterdodge(jitter.width = 0.0, dodge.width =
                                               0.75),
             size = 1.5) +
  stat_summary(
    fun = mean,
    position = position_dodge(0.75),
    geom = "point",
    size = 2.5,
    pch = 23
  ) +
  scale_color_viridis_c(option = "C", name = TeX(T0string)) +
  scale_fill_viridis_d(option = "mako", name = TeX("SiO$_2$ (wt%)")) +
  labs(
    y = TeX("Densification rate (kg/m$^3$/Myr)"),
    x = TeX("Damköhler number")
  ) +
  #annotate("text",x=0.1,y=1.28e3,label="reaction limited",colour="#666666") +
  #annotate("segment",x=0.003,xend=100,y=1.08e3,yend=1.08e3,colour='#999999',size=1, arrow=arrow(type="closed",length=unit(0.01,"npc"),ends="both")) +
  #annotate("text",x=6e5,y=1.28e3,label="density limited",colour="#666666") +
  #annotate("segment",x=1e4,xend=3e6,y=1.08e3,yend=1.08e3,colour='#999999',size=1, arrow=arrow(type="closed",length=unit(0.01,"npc"),ends="both")) +
  theme(legend.position = "bottom",
        legend.title = element_text(vjust = 0.75))

ggsave(
  'plots/densification_rate.pdf',
  device = cairo_pdf,
  width = 6.5,
  height = 6
)
ggsave('plots/densification_rate.png',
       width = 6.5,
       height = 6)



dat %>%
  ggplot(aes(x = r0 * 3.154e7 * 1e6 , y = 2780 / (r0 * 6000 * 3.154e7))) +
  #geom_abline(slope=0,intercept=1) +
  #geom_abline(slope=0,intercept=log10(50)) +
  geom_hline(yintercept = 3e6,
             color = "#555555",
             linetype = "63") +
  geom_hline(yintercept = 3.23e6,
             size = 26,
             color = "#ff333366") +
  geom_vline(
    xintercept = 1.5e-10 * 1e6 / 1e3 * 100 * 100,
    size = 20,
    colour = "#0066aa66"
  ) +
  geom_vline(
    xintercept = 1.5e-9 * 1e6 / 1e3 * 100 * 100,
    size = 20,
    colour = "#33333333"
  ) +
  geom_vline(
    xintercept = 2.3e-4 * 1e6 / 1e3 * 100 * 100,
    size = 20,
    colour = "#00aa3366"
  ) +
  geom_vline(xintercept = 0.15, linetype = "63") +
  labs(
    x = TeX("Eff. rate constant $\\r_0$ (kg/m$^2$/Myr)"),
    y = TeX("Reactive timescale $\\tau$ (yr)")
  ) +
  geom_smooth(method = lm,
              se = FALSE,
              fullrange = TRUE) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4),
                labels = math_format(.x)) +
  scale_y_log10(
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12),
    labels =  trans_format('log10', math_format(10 ^ .x))
  ) +
  theme_mcm()


selected_settings <-
  c(
    'hot-1',
    'hot-2',
    'hot-3',
    'transitional-1',
    'transitional-2',
    'transitional-3',
    'cold-1',
    'cold-2',
    'cold-3'
  )

for (i in seq_along(selected_settings)) {
  selected_setting <- selected_settings[i]
  print(selected_setting)
  
  drho1 <- dat %>%
    filter(
      composition == 'mackwell_1998_maryland_diabase_norm' &
        setting == selected_setting & Da == 1000
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho2 <- dat %>%
    filter(
      composition == 'mackwell_1998_maryland_diabase_norm' &
        setting == selected_setting & Da == 100
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho3 <- dat %>%
    filter(
      composition == 'mackwell_1998_maryland_diabase_norm' &
        setting == selected_setting & Da == 10
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho4 <- dat %>%
    filter(
      composition == 'mackwell_1998_maryland_diabase_norm' &
        setting == selected_setting & Da == 1
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  
  drho5 <- dat %>%
    filter(
      composition == 'hacker_2015_md_xenolith_norm' &
        setting == selected_setting & Da == 1000
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho6 <- dat %>%
    filter(
      composition == 'hacker_2015_md_xenolith_norm' &
        setting == selected_setting & Da == 100
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho7 <- dat %>%
    filter(composition == 'hacker_2015_md_xenolith_norm' &
             setting == selected_setting & Da == 10) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho8 <- dat %>%
    filter(composition == 'hacker_2015_md_xenolith_norm' &
             setting == selected_setting & Da == 1) %>%
    select(effective_delta_rho) %>% pull * 1000
  
  drho9 <- dat %>%
    filter(
      composition == 'sammon_2021_lower_crust_norm' &
        setting == selected_setting & Da == 1000
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho10 <- dat %>%
    filter(
      composition == 'sammon_2021_lower_crust_norm' &
        setting == selected_setting & Da == 100
    ) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho11 <- dat %>%
    filter(composition == 'sammon_2021_lower_crust_norm' &
             setting == selected_setting & Da == 10) %>%
    select(effective_delta_rho) %>% pull * 1000
  drho12 <- dat %>%
    filter(composition == 'sammon_2021_lower_crust_norm' &
             setting == selected_setting & Da == 1) %>%
    select(effective_delta_rho) %>% pull * 1000
  
  print(drho1)
  print(drho2)
  print(drho3)
  print(drho4)
  print(drho5)
  print(drho6)
  print(drho7)
  print(drho8)
  print(drho9)
  print(drho10)
  print(drho11)
  print(drho12)
  # Jull & Kelemen using Hirth 2003 dry olivine disl.
  
  Tlow = dat %>% filter(setting == selected_setting) %>% select(T0) %>% unique %>% pull
  Thigh = dat %>% filter(setting == selected_setting) %>% select(T1) %>% unique %>% pull
  Temp = Tlow / 4 + 3 * Thigh / 4
  Temp = Thigh
  tb <- function(h_km, drho) {
    if (is.nan(drho)) {
      return(NaN)
    }
    
    # h in meters
    h = h_km * 1000
    # drho in kg/m3
    
    R = 8.3145 # J/mol/K
    g = 9.81 # m/s2
    n = 3.5
    
    Q = 535e3 # J/mol, wet olivine
    
    #Q = 530e3
    #A = 1.1e5 # MPa^(-n) s^(-1)
    
    # need to better constrain these
    z0prime = 0.01 # initial displacement in layer thicknesses
    zprime = 1 # downward displacement in layer thicknesses
    qprime = 0.3 # growth rate
    
    # Newtonian
    n = 1
    A_Mpa = 1.5e9 # Mpa^-n s^-1
    #A_Mpa = exp(15.4)
    Q = 535e3
    A = A_Mpa * 1e6 ^ (-n) # Pa^-n s^-1
    B = A ^ (-1 / n) * exp(Q / (n * R * Temp)) # Pa s = kg/m/s
    Tb = (B / (2 * drho * g * h)) ^ n # s
    tbn = (-1 / qprime) * log(z0prime / zprime) * Tb # seconds
    tbn = tbn / 3.154e7 # yrs
    
    # non-Newtonian
    n = 3.5
    Cprime = 20 # fig 15
    A_Mpa = exp(15.4) # MPa^-n s^-1, wet olivine
    A = A_Mpa * 1e6 ^ (-n) # Pa^-n s^-1
    Q = 515e3
    B = A ^ (-1 / n) * exp(Q / (n * R * Temp)) # Pa s = kg/m/s
    Tb = (B / (2 * drho * g * h)) ^ n # s
    tbnn = (n / Cprime) ^ n * (z0prime ^ (1 - n)) / (n - 1) * Tb # seconds
    tbnn = tbnn / 3.154e7 # yr
    
    strain_rate = 1e-15
    eprime = strain_rate * Tb
    
    eprime0 = 1e-10
    dt = 2e5 - 2e7
    de = 10e-6 - 10e-10
    dtde = dt / de
    dtde = -0.5
    tbnn = (tbnn) * (eprime / eprime0) ^ -(eprime / tbnn * dtde / tbnn)
    
    return(tbnn)
  }
  
  tyr <- 10 ^ (seq(4, 10, .01)) # 10,000 to 100,000,000
  
  rate <- 1e-3 # m/yr
  hs <- tyr * rate / 1000
  
  tryCatch(
    expr = {
      theplot <- ggplot() +
        geom_ribbon(
          aes(
            ymin = tb(hs, drho12),
            ymax = tb(hs, drho9),
            x = hs,
            fill = 'A'
          ),
          alpha = 0.1,
          show.legend = FALSE
        ) +
        geom_ribbon(
          aes(
            ymin = tb(hs, drho8),
            ymax = tb(hs, drho5),
            x = hs,
            fill = 'B'
          ),
          alpha = 0.1,
          show.legend = FALSE
        ) +
        
        geom_ribbon(
          aes(
            ymin = tb(hs, drho4),
            ymax = tb(hs, drho1),
            x = hs,
            fill = 'C'
          ),
          alpha = 0.1,
          show.legend = FALSE
        ) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho1),
          colour = 'C',
          linetype = 'A',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho2),
          colour = 'C',
          linetype = 'B',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho3),
          colour = 'C',
          linetype = 'C',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho4),
          colour = 'C',
          linetype = 'D',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho5),
          colour = 'B',
          linetype = 'A',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho6),
          colour = 'B',
          linetype = 'B',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho7),
          colour = 'B',
          linetype = 'C',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho8),
          colour = 'B',
          linetype = 'D',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho9),
          colour = 'A',
          linetype = 'A',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho10),
          colour = 'A',
          linetype = 'B',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho11),
          colour = 'A',
          linetype = 'C',
        ), size = 0.5) +
        geom_line(aes(
          x = hs,
          y = tb(hs, drho12),
          colour = 'A',
          linetype = 'D',
        ), size = 0.5) +
        coord_cartesian(xlim = c(0, 40),
                        ylim = c(1e5, 5e7)) +
        scale_y_continuous(
          expand = c(0, 0),
          trans = "log10",
          breaks = c(1e5, 1e6, 1e7, 1e8, 1e9),
          labels = trans_format('log10', math_format(10 ^ .x))
        ) +
        annotation_logticks(
          sides = "l",
          size = 0.25,
          short = unit(0.1, "cm"),
          mid = unit(0.2, "cm"),
          long = unit(0.3, "cm")
        ) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5)) +
        labs(x = "Layer thickness (km)", y = "Time (yr)") +
        scale_linetype_manual(
          values = c("solid", "solid", "solid", "solid"),
          labels = c("100", "10", "3", "1"),
          name = "Da",
        ) +
        guides(linetype = guide_legend(reverse = T)) +
        scale_colour_discrete(name = "Composition") +
        geom_ribbon(
          aes(
            xmin = hs * 0.7,
            xmax = hs * 1.5,
            y = tyr
          ),
          fill = "#dadada",
          alpha = 0.4
        ) +
        geom_line(aes(x = hs, y = tyr),
                  colour = '#333333',
                  size = 1) +
        theme_mcm() +
        theme(legend.position = 'none')
      
      
      
      ggsave(
        paste('plots/thickness_', selected_setting, '.pdf', sep = ""),
        device = cairo_pdf,
        width = 3.75,
        height = 3.5
      )
      ggsave(
        paste('plots/thickness_', selected_setting, '.png', sep = ""),
        width = 3.75,
        height = 3.5
      )
      print(theplot)
    },
    error = function(cond) {
      message(conditionMessage(cond))
      # Choose a return value in case of error
      NA
    }
  )
  
}
