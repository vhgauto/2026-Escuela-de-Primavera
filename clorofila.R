# Preprocesamiento

library(terra)
library(tidyterra)
library(tidymodels)
library(tidyverse)
tidymodels_prefer()

l <- list.files("tp_clorofila/raster_salida/", full.names = TRUE)
r_lista <- map(l, rast)

p <- vect("tp_clorofila/vectores/puntos.geojson")

r2 <- r_lista[[2]]
r2_fecha <- r2 |>
  sources() |>
  basename() |>
  ymd()

ggplot() +
  geom_spatraster_rgb(
    data = stretch(r2),
    r = 4,
    g = 3,
    b = 2,
    stretch = "lin"
  ) +
  geom_spatvector(data = p, shape = 21, fill = "red", color = "gold") +
  geom_spatvector_label(
    data = p,
    aes(label = sitio),
    size = 2,
    vjust = -.25,
    border.color = NA,
    fill = "white"
  ) +
  labs(title = r2_fecha, x = NULL, y = NULL) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 8) +
  theme_sub_plot(background = element_blank()) +
  theme_sub_panel(background = element_blank())

p <- p[p$sitio != "GAR"]

f_reflect <- function(i) {
  terra::extract(i, p, bind = TRUE) |>
    as_tibble() |>
    mutate(fecha = ymd(basename(sources(i))), .before = 1)
}

reflect_tbl <- map(r_lista, f_reflect) |>
  list_rbind() |>
  select(fecha, sitio, starts_with("B"))

cla_tbl <- read_csv("tp_clorofila/datos/clorofila.csv")

reflect_cla <- inner_join(reflect_tbl, cla_tbl, by = join_by(fecha, sitio)) |>
  select(cla, starts_with("B"))

# Correlaciones

d <- reflect_cla |>
  mutate(
    ndvi = (B5 - B4) / (B5 + B4),
    B5_B4 = B5 / B4
  )

pivot_longer(
  d,
  cols = -cla,
  values_to = "reflect",
  names_to = "banda"
) |>
  nest(.by = banda) |>
  mutate(r = map_dbl(data, ~ cor(.x$reflect, .x$cla))) |>
  mutate(pvalor = map(data, ~ cor.test(.x$reflect, .x$cla))) |>
  mutate(pvalor = map_dbl(pvalor, ~ tidy(.x)$p.value)) |>
  select(banda, r, pvalor) |>
  arrange(desc(r)) |>
  mutate(es_signif = pvalor < .05)

# División de Datos

set.seed(2000)
datos_split <- initial_split(d, strata = cla, prop = 3 / 4)
datos_train <- training(datos_split)
datos_test <- testing(datos_split)
datos_split

# Re-muestreo

set.seed(2001)
datos_bootstrap <- bootstraps(datos_train, strata = cla, times = 25)
datos_bootstrap

# Construcción de Modelos

## Mecanismos de los modelos

lin_spec <- linear_reg(mode = "regression") |>
  set_engine("lm")

rf_spec <- rand_forest(trees = tune()) |>
  set_engine("ranger") |>
  set_mode("regression")

## Recetas

receta_ndvi <- recipe(cla ~ ndvi, data = d)
receta_b5 <- recipe(cla ~ B5, data = d)
receta_b5_b4 <- recipe(cla ~ B5_B4, data = d)

receta_bandas <- recipe(cla ~ B1 + B2 + B3 + B4 + B5 + B6 + B7, data = d)
receta_select <- recipe(cla ~ B5 + B3, data = d)

## Workflow

wf <- workflow_set(
  preproc = list(
    ndvi = receta_ndvi,
    b5 = receta_b5,
    b5_b4 = receta_b5,
    bandas = receta_bandas,
    select = receta_select
  ),
  models = list(
    lin = lin_spec,
    lin = lin_spec,
    lin = lin_spec,
    rf = rf_spec,
    rf = rf_spec
  ),
  cross = FALSE
)
wf

# Entrenamiento

race_ctrl <- finetune::control_race(
  save_pred = TRUE,
  parallel_over = "everything",
  save_workflow = TRUE
)

res <- wf |>
  workflow_map(
    fn = "tune_race_anova",
    seed = 2002,
    resamples = datos_bootstrap,
    grid = 50,
    control = race_ctrl
  )

# Selección del modelo

f_entrenamiento <- function(metrica) {
  autoplot(
    res,
    rank_metric = metrica,
    metric = metrica,
    select_best = TRUE
  ) +
    geom_text(
      aes(y = mean, label = wflow_id),
      angle = 90,
      vjust = -.3,
      size = 5
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(aspect.ratio = 1, legend.position = "none") +
    theme_sub_panel(background = element_blank()) +
    theme_sub_plot(background = element_blank())
}
f_entrenamiento("rmse")
f_entrenamiento("rsq")

res |>
  rank_results() |>
  pivot_wider(
    names_from = .metric,
    values_from = mean,
    id_cols = c(wflow_id, .config, rank, model)
  ) |>
  slice_min(order_by = rmse, n = 1, with_ties = FALSE) |>
  select(wflow_id, rank, rmse, rsq)

mejor_res <- res |>
  extract_workflow_set_result("ndvi_lin") |>
  select_best(metric = "rmse")

# Validación

val_res <- res |>
  extract_workflow("ndvi_lin") |>
  finalize_workflow(mejor_res) |>
  last_fit(split = datos_split)

collect_predictions(val_res) |>
  ggplot(aes(x = cla, y = .pred)) +
  geom_abline(color = "gray50", lty = 2) +
  geom_point(size = 4, alpha = .5, color = "dodgerblue") +
  annotate(geom = "text", x = I(.9), y = I(.94), label = "y=x", angle = 45) +
  coord_obs_pred(expand = TRUE) +
  labs(x = "cla observado (μg/L)", y = "cla estimado (μg/L)") +
  theme_bw() +
  theme_sub_panel(
    grid.major = element_line(linewidth = .1),
    background = element_blank()
  ) +
  theme_sub_plot(background = element_blank())

collect_metrics(val_res) |>
  pivot_wider(names_from = .metric, values_from = .estimate) |>
  select(rmse, rsq)

collect_predictions(val_res) |>
  lm(.pred ~ cla, data = _) |>
  confint(level = .95)

const <- extract_fit_engine(val_res) |>
  tidy() |>
  pull(estimate) |>
  round(2)

cat("cla =", const[1], "+", const[2], "* NDVI")

# Mapa

mndwi <- (r2$B3 - r2$B6) / (r2$B3 + r2$B6)

m <- thresh(mndwi, method = "otsu")
m[isTRUE(m)] <- 1
m[isFALSE(m)] <- NA

agua <- r2 * m
agua_ndvi <- (agua$B5 - agua$B4) / (agua$B5 + agua$B4)
terra::set.names(agua_ndvi, "ndvi")

plot(m, axes = FALSE, col = "dodgerblue", legend = FALSE, box = TRUE)
north(xy = "topleft", type = 1)

cla_pred <- terra::predict(agua_ndvi, extract_fit_engine(val_res))
cla_pred <- clamp(cla_pred, 0, 80)

ggplot() +
  geom_spatraster(data = cla_pred) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    width = unit(20, "pt"),
    height = unit(30, "pt")
  ) +
  ggspatial::annotation_scale(location = "br", width_hint = .3) +
  scale_fill_whitebox_c(palette = "muted", breaks = scales::breaks_width(10)) +
  labs(title = r2_fecha, fill = "Clorofila-a\n(μg/L)") +
  theme_minimal() +
  theme_sub_legend(key.height = unit(45, "pt")) +
  theme_sub_panel(
    background = element_rect(fill = "grey95"),
    grid.major = element_line(linetype = 2, color = "grey80", linewidth = .3)
  ) +
  theme_sub_plot(background = element_blank()) +
  theme_sub_axis_left(text = element_text(angle = 90, hjust = .5))

# Explicabilidad

rf_mejor_res <- res |>
  extract_workflow_set_result("select_rf") |>
  select_best(metric = "rmse")

rf_val_res <- res |>
  extract_workflow("select_rf") |>
  finalize_workflow(rf_mejor_res) |>
  last_fit(split = datos_split)

vip_train <- datos_train |>
  select(B3, B5)

explainer <- DALEXtra::explain_tidymodels(
  model = extract_workflow(rf_val_res),
  data = vip_train,
  y = datos_train$cla,
  label = "RF",
  verbose = FALSE
)

set.seed(2026)
vip <- DALEX::model_parts(
  explainer = explainer,
  loss_function = DALEX::loss_root_mean_square
)

vip_df <- as_tibble(vip) |>
  filter(variable != "_full_model_") |>
  filter(variable != "_baseline_") |>
  mutate(variable = fct_reorder(variable, dropout_loss))

ggplot(vip_df, aes(dropout_loss, variable)) +
  geom_boxplot() +
  labs(y = NULL, x = "RMSE (μg/L)") +
  theme_bw() +
  theme_sub_panel(background = element_blank()) +
  theme_sub_plot(background = element_blank())
