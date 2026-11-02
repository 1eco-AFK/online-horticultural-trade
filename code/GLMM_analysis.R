# GLMM fitting, model selection, diagnostics and likelihood-ratio tests.
# Run from the repository root. This script fits models from the trade CSV.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lme4)
  library(MuMIn)
  library(performance)
})
dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
options(na.action = "na.fail")

# 1. Aggregate records and construct the original orthogonal contrasts.
df <- read_csv("data/Step4_Refugia_Full_Filtered_NotVeg.csv", show_col_types = FALSE)
glimpse(df)
df_agg <- df %>%
  group_by(shop_id, wcg_lat, wcg_lng, year_label, invasion_label, propagule_type) %>%
  summarise(
    refugia_count = sum(is_refugia),
    total_count = n(),
    ratio = refugia_count / total_count,
    not.refugia_count = total_count - refugia_count,
    .groups = 'drop'
  ) %>%
  mutate(
    shop_id = as.factor(shop_id),
    year_label = as.factor(year_label),
    invasion_label = as.factor(invasion_label),
    propagule_type = as.factor(propagule_type)
  )
cat("Aggregated Data Summary:\n")
print(head(df_agg))
write_csv(df_agg, "results/Step4_df_agg.csv")
mat_comp <- cbind(
  c(1/3, 1/3, 1/3),
  c(-1, 0.5, 0.5),
  c(0, -1, 1)
)
mymat_comp <- solve(t(mat_comp))
contrast_comp <- as.data.frame(mymat_comp[, 2:3])
contrast_comp$invasion_label <- c("nonnat", "nat", "inv")
names(contrast_comp)[1:2] <- c("C1_NatStatus", "C2_Invasiveness")
df_agg <- df_agg %>% left_join(contrast_comp, by = "invasion_label")
glimpse(df_agg)
df_agg$wcg_lat_scaled <- scale(df_agg$wcg_lat)
df_agg$wcg_lng_scaled <- scale(df_agg$wcg_lng)
strict_control <- glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5))
model_interaction_full <- glmer(cbind(refugia_count, not.refugia_count) ~
                                  wcg_lat_scaled +
                                  wcg_lng_scaled +
                                  C1_NatStatus +
                                  C2_Invasiveness +
                                  propagule_type +
                                  year_label +
                                  wcg_lat_scaled:C1_NatStatus +
                                  wcg_lat_scaled:C2_Invasiveness +
                                  wcg_lat_scaled:propagule_type +
                                  wcg_lat_scaled:year_label +
                                  wcg_lng_scaled:C1_NatStatus +
                                  wcg_lng_scaled:C2_Invasiveness +
                                  wcg_lng_scaled:propagule_type +
                                  wcg_lng_scaled:year_label +
                                  (1 | shop_id),
                                data = df_agg,
                                family = binomial,
                                control = strict_control)

stopifnot(nrow(df) == 21351L, nrow(df_agg) == 3333L)
write_csv(mutate(df_agg, across(ends_with("_scaled"), as.numeric)),
          "results/Step4_model_data.csv")

# 2. Model selection.
n_cores <- 6L
worker_cluster <- parallel::makeCluster(n_cores)
parallel::clusterExport(worker_cluster, c("df_agg", "strict_control"), envir = .GlobalEnv)
parallel::clusterEvalQ(worker_cluster, {library(lme4); options(na.action = "na.fail")})
message("Fitting candidate models with ", n_cores, " workers.")
dd <- MuMIn::dredge(model_interaction_full, cluster = worker_cluster)
parallel::stopCluster(worker_cluster)
best_model <- MuMIn::get.models(dd, subset = 1)[[1]]
saveRDS(model_interaction_full, "models/Step4_Full_Model_Results.rds")
saveRDS(dd, "models/Step4_Dredge_Results.rds")
saveRDS(best_model, "models/Step4_best_model_Results.rds")
write_csv(as.data.frame(dd), "results/Step4_Model_Selection.csv")

# 3. Report estimates and diagnostics for the full and selected models.
coefficient_table <- as.data.frame(coef(summary(best_model)))
coefficient_table$term <- rownames(coefficient_table)
write_csv(coefficient_table, "results/Step4_Model_Coefficients.csv")
sink("results/Step4_Model_Diagnostics.txt")
print(summary(model_interaction_full))
print(model_interaction_full@optinfo$conv)
print(performance::check_overdispersion(model_interaction_full))
print(performance::check_collinearity(model_interaction_full))
print(car::vif(model_interaction_full))
print(summary(best_model))
print(best_model@optinfo$conv)
cat("Singular full model:", lme4::isSingular(model_interaction_full), "\n")
cat("Singular selected model:", lme4::isSingular(best_model), "\n")
print(performance::check_overdispersion(best_model))
print(performance::check_collinearity(best_model))
sink()
write_csv(as.data.frame(performance::model_performance(best_model)),
          "results/Step4_Model_Performance_Metrics.csv")

# 4. Interaction tests.
m_no_int1 <- update(best_model, . ~ . - C1_NatStatus:wcg_lat_scaled)
lrt_int1 <- anova(m_no_int1, best_model, test = "Chisq")

m_no_int2 <- update(best_model, . ~ . - C1_NatStatus:wcg_lng_scaled)
lrt_int2 <- anova(m_no_int2, best_model, test = "Chisq")

m_no_int3 <- update(best_model, . ~ . - C2_Invasiveness:wcg_lat_scaled)
lrt_int3 <- anova(m_no_int3, best_model, test = "Chisq")

m_no_int4 <- update(best_model, . ~ . - propagule_type:wcg_lat_scaled)
lrt_int4 <- anova(m_no_int4, best_model, test = "Chisq")

# 5. Main-effect tests.
m_additive <- update(best_model, . ~ . - C1_NatStatus:wcg_lat_scaled
                                       - C1_NatStatus:wcg_lng_scaled
                                       - C2_Invasiveness:wcg_lat_scaled
                                       - propagule_type:wcg_lat_scaled)

m_no_main1 <- update(m_additive, . ~ . - C1_NatStatus)
lrt_main1 <- anova(m_no_main1, m_additive, test = "Chisq")

m_no_main2 <- update(m_additive, . ~ . - C2_Invasiveness)
lrt_main2 <- anova(m_no_main2, m_additive, test = "Chisq")

m_no_main3 <- update(m_additive, . ~ . - propagule_type)
lrt_main3 <- anova(m_no_main3, m_additive, test = "Chisq")

m_no_main4 <- update(m_additive, . ~ . - wcg_lat_scaled)
lrt_main4 <- anova(m_no_main4, m_additive, test = "Chisq")

m_no_main5 <- update(m_additive, . ~ . - wcg_lng_scaled)
lrt_main5 <- anova(m_no_main5, m_additive, test = "Chisq")

m_no_main6 <- update(m_additive, . ~ . - year_label)
lrt_main6 <- anova(m_no_main6, m_additive, test = "Chisq")

# 6. Save the test results.
lrt_table <- data.frame(
  term = c("C1_NatStatus:wcg_lat_scaled", "C1_NatStatus:wcg_lng_scaled",
           "C2_Invasiveness:wcg_lat_scaled", "propagule_type:wcg_lat_scaled",
           "C1_NatStatus", "C2_Invasiveness", "propagule_type",
           "wcg_lat_scaled", "wcg_lng_scaled", "year_label"),
  reference_model = c(rep("selected", 4), rep("additive", 6)),
  AIC_full = c(rep(AIC(best_model), 4), rep(AIC(m_additive), 6)),
  AIC_reduced = c(AIC(m_no_int1), AIC(m_no_int2), AIC(m_no_int3), AIC(m_no_int4),
                  AIC(m_no_main1), AIC(m_no_main2), AIC(m_no_main3),
                  AIC(m_no_main4), AIC(m_no_main5), AIC(m_no_main6)),
  Chisq = c(lrt_int1$Chisq[2], lrt_int2$Chisq[2], lrt_int3$Chisq[2], lrt_int4$Chisq[2],
            lrt_main1$Chisq[2], lrt_main2$Chisq[2], lrt_main3$Chisq[2],
            lrt_main4$Chisq[2], lrt_main5$Chisq[2], lrt_main6$Chisq[2]),
  Df = c(lrt_int1$Df[2], lrt_int2$Df[2], lrt_int3$Df[2], lrt_int4$Df[2],
         lrt_main1$Df[2], lrt_main2$Df[2], lrt_main3$Df[2],
         lrt_main4$Df[2], lrt_main5$Df[2], lrt_main6$Df[2]),
  P_value = c(lrt_int1$`Pr(>Chisq)`[2], lrt_int2$`Pr(>Chisq)`[2],
              lrt_int3$`Pr(>Chisq)`[2], lrt_int4$`Pr(>Chisq)`[2],
              lrt_main1$`Pr(>Chisq)`[2], lrt_main2$`Pr(>Chisq)`[2],
              lrt_main3$`Pr(>Chisq)`[2], lrt_main4$`Pr(>Chisq)`[2],
              lrt_main5$`Pr(>Chisq)`[2], lrt_main6$`Pr(>Chisq)`[2])
)
write_csv(lrt_table, "results/Step4_Model_LRT.csv")
capture.output(sessionInfo(), file = "results/GLMM_sessionInfo.txt")
message("GLMM workflow complete. Run code/Figure_2.R and code/Figure_3.R next.")

