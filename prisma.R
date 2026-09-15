# Obtención de Datos

library(terra)
library(tidyverse)
library(tidyterra)

r <- rast("tp_prisma/salida/20220311.tif")
nlyr(r)

p <- vect("tp_prisma/vectores/puntos.geojson")

r_fecha <- sources(r) |>
  basename() |>
  ymd()

ggplot() +
  geom_spatraster_rgb(
    data = r,
    r = 34,
    g = 20,
    b = 10,
    stretch = "lin"
  ) +
  geom_spatvector(data = p, shape = 21, fill = "red", color = "gold") +
  geom_spatvector_label(
    data = p,
    aes(label = puntos),
    size = 2,
    vjust = -.25,
    border.color = NA,
    fill = "white"
  ) +
  labs(title = r_fecha, x = NULL, y = NULL) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 8) +
  theme_sub_plot(background = element_blank()) +
  theme_sub_panel(background = element_blank())

# Firmas espectrales

wvl <- read_delim(
  "tp_prisma/salida/PRS_L2D_STD_20220311142257_20220311142301_0001_HCO_FULL.wvl",
  delim = " ",
  show_col_types = FALSE
) |>
  rename(banda = band) |>
  select(wl, banda)

reflect_tbl <- terra::extract(r, p, bind = TRUE) |>
  as_tibble() |>
  pivot_longer(
    cols = starts_with("PRS_L2D_STD_"),
    names_to = "banda",
    values_to = "reflect"
  ) |>
  mutate(banda = str_extract(banda, "\\d+$")) |>
  mutate(banda = as.numeric(banda)) |>
  inner_join(wvl, by = join_by(banda)) |>
  select(puntos, banda, reflect, wl)

ggplot(reflect_tbl, aes(wl, reflect, group = puntos, color = puntos)) +
  geom_line(linewidth = .5) +
  scale_x_continuous(breaks = scales::breaks_width(200)) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Long. de onda (nm)", y = "Reflect.", color = NULL) +
  theme_bw() +
  theme_sub_axis(text = element_text(color = "black")) +
  theme_sub_panel(
    grid.minor = element_blank(),
    grid.major = element_line(linewidth = .1, color = "grey"),
    background = element_blank()
  ) +
  theme_sub_strip(text = element_text(color = "black", face = "bold")) +
  theme_sub_plot(background = element_blank(), margin = margin(r = 10)) +
  theme_sub_legend(
    position = "top",
    key.spacing.x = unit(20, "pt"),
    background = element_blank(),
    text = element_text(margin = margin())
  )

# Spectral Angle Mapper (SAM)

banda_rango <- wvl |>
  filter(between(wl, 400, 1800)) |>
  pull(banda) |>
  range()
banda_rango

r_rango <- r[[banda_rango[1]:banda_rango[2]]]
reflect_vnir <- terra::extract(r_rango, p)
rownames(reflect_vnir) <- p$puntos

lago_sam_ang <- RStoolbox::sam(r_rango, reflect_vnir, angles = TRUE)

plot(
  lago_sam_ang,
  nc = 2,
  nr = 2,
  axes = FALSE,
  box = TRUE,
  col = rev(viridis::magma(500))
)

lago_sam_clas <- RStoolbox::sam(r_rango, reflect_vnir, angles = FALSE)

ggplot() +
  geom_spatraster(data = lago_sam_clas) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    width = unit(20, "pt"),
    height = unit(30, "pt")
  ) +
  ggspatial::annotation_scale(location = "bl", width_hint = .3) +
  scale_fill_princess_c(
    palette = "maori",
    breaks = 1:4,
    labels = p$puntos,
    direction = -1,
    guide = "legend"
  ) +
  labs(fill = NULL) +
  coord_sf(expand = FALSE) +
  labs(title = r_fecha) +
  theme_minimal(base_size = 8) +
  theme_sub_legend(
    position = "top",
    key.spacing.x = unit(15, "pt"),
    text = element_text(margin = margin(l = 3))
  ) +
  theme_sub_panel(
    background = element_rect(fill = "grey95"),
    grid.major = element_line(linetype = 2, color = "grey80", linewidth = .3)
  ) +
  theme_sub_plot(background = element_blank()) +
  theme_sub_axis_left(text = element_text(angle = 90, hjust = .5))

# Principal Component Analysis (PCA)

mndwi <- (r[[34]] - r[[179]]) / (r[[34]] + r[[179]])
m <- thresh(mndwi)
m[isTRUE(m)] <- 1
m[isFALSE(m)] <- NA
m <- fillHoles(m)

plot(m, col = "dodgerblue", axes = FALSE, box = TRUE, legend = FALSE)
north(xy = "topleft")

agua_pca <- m * r_rango
r_pca <- RStoolbox::rasterPCA(agua_pca, nComp = 3, spca = TRUE)

plot(
  r_pca$map,
  col = viridis::magma(500),
  axes = FALSE,
  box = TRUE,
  nc = 3,
  nr = 1
)

plotRGB(stretch(r_pca$map), stretch = "lin", colNA = "transparent")
north(xy = "topleft")

acumulado_pca <- factoextra::get_eigenvalue(r_pca$model) |>
  as_tibble(rownames = "dim") |>
  slice_head(n = 3) |>
  select(dim, var = variance.percent, acumulado = cumulative.variance.percent)
acumulado_pca

acumulado_pca |>
  mutate(acumulado_label = paste0(round(acumulado, 1), "%")) |>
  ggplot(aes(x = dim, fill = dim)) +
  geom_col(aes(y = var)) +
  geom_line(aes(y = acumulado, group = 1)) +
  geom_point(
    aes(y = acumulado, fill = dim),
    shape = 21,
    size = 3,
    stroke = 1,
    color = "white"
  ) +
  geom_text(aes(y = acumulado, label = acumulado_label), vjust = -1) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  scale_fill_brewer(palette = "Dark2", guide = guide_none()) +
  coord_cartesian(clip = "off", ylim = c(0, 100)) +
  labs(x = NULL, y = "Aporte (%)") +
  theme_bw() +
  theme_sub_axis(text = element_text(color = "black")) +
  theme_sub_plot(background = element_blank()) +
  theme_sub_panel(background = element_blank())
