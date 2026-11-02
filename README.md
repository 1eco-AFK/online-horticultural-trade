# Plant trade and climatic refugia

`code/`: R scripts; `data/`: inputs; `models/`: saved models and selection results; `results/`: Figures 1-3 and S1-S2 (PNG/PDF).

Use R 4.4.1 with Cairo support. Set the working directory to this folder.

```r
install.packages(c("tidyverse", "sf", "terra", "ggsci", "ggpubr", "patchwork",
                   "lme4", "MuMIn", "performance", "car", "ggeffects", "emmeans", "ragg"))

# Draw all five figures using the supplied models and S2 cache.
source("code/Figure_1.R")
source("code/Figure_2.R")
source("code/Figure_3.R")
source("code/Figure_S1.R")
source("code/Figure_S2.R")
```

To refit models from the trade data, run `source("code/GLMM_analysis.R")` before drawing Figures 2-3. To rebuild S2 data from source rasters, run `source("code/Figure_S2_data.R")` before `source("code/Figure_S2.R")`.
