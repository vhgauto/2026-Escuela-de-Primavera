# Escuela de Primavera 2026

**Remediación, monitoreo y gestión de cuencas antropizadas: Lecciones aprendidas de la teledetección**

Más detalles en el [sitio del evento](https://ig.conae.unc.edu.ar/escuela-de-primavera-2026/).

---

El presente repositorio contiene los datos y scripts necesarios para las actividades prácticas de modelización de clorofila-a (lunes 21/9) y procesamiento de una imagen hiperespectral (martes 22/9).

Para facilitar las actividades durante la clase es conveniente clonar este repositorio y trabajar de manera local.

El lenguaje de programación utilizado en ambas actividades prácticas es `R`. Se recomienda tenerlo [instalado](https://cloud.r-project.org/) y contar con un editor, como ser [RStudio](https://docs.posit.co/ide/user/#rstudio-ide-oss-downloads) o [Positron](https://positron.posit.co/download.html).

La instalación de los paquetes está organizada por [`renv`](https://rstudio.github.io/renv/articles/renv.html). Pueden instalar todos los paquetes necesarios para ambas actividades utilizando el siguiente comando.

```r
# install.packages("renv")
renv::restore()
```

Hay disponibles dos scripts (`clorofila.R` y `prisma.R`) para quienes deseen seguir la clase, ejecutando paso a paso los comandos. Asimismo, los scripts también se encuentran disponibles como notebooks en Google Colab ([clorofila-a](https://colab.research.google.com/drive/1yS79LNB4eCzfhuDCTfH5m-2peRXGGdY8) y [PRISMA](https://colab.research.google.com/drive/1d8UzcB0-F6MLvQQDFCriV_T8c5nOBFFS)).

## Integrantes

- Dra. Anabella Ferral
- Dra. Rocío Guido
- Mgtr. Víctor Gauto
