## ===========================================
###Spatial TME risk score related analysis###
## ===========================================
library(survminer)
library(survival)
library(glmnet)
library(tidyverse)
library(tidyr)
library(dplyr)
library(forestplot)
library(tableone)
library(plotmo)
library(pROC)
library(car)
setwd('D:/BaiduSyncdisk/iCCA/cell_frequency')
spe_ICC <- readRDS('spe_ICC.rds')
clinical_155 <- read.csv('clinical_155.csv')
## Preprocessing: Compute cell frequencies
# Extract image ID and cell annotation from spe_ICC
a <- spe_ICC@colData$ImageNb
b <- spe_ICC@colData$detailed_anno
c <- as.data.frame(cbind(a, b))

# Count the number of each cell type per image
result <- c %>%
  group_by(a, b) %>%
  summarise(count = n()) %>%
  spread(b, count, fill = 0)

print(result)

# Convert to row-based form and calculate percentage frequency of each cell type
result <- column_to_rownames(result, var = 'a')
row_sums <- rowSums(result)
ratio <- sweep(result, 1, row_sums, "/") * 100
ratio <- rownames_to_column(ratio, var = 'a')

# Merge with clinical information (assumes roi_id in column 1)
clinical <- clinical_155[, c(1, 15, 16)]
result_2 <- merge(ratio, clinical, by.x = 'a', by.y = 'roi_id')

# write.csv(result_2, 'TME_cell_frequency_cli.CSV')

## Univariate Cox regression: filter significant predictors
TME <- read.csv('TME_cell_frequency_cli.csv', check.names = FALSE, row.names = 1)
coxPfilter <- 0.05
outTab <- data.frame()
sigGenes <- c("OS", "OS_event")

for (i in colnames(TME[, 1:37])) {
  cox <- coxph(Surv(OS, OS_event) ~ TME[, i], data = TME)
  coxSummary <- summary(cox)
  coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
  
  # Keep significant variables
  if (coxP < coxPfilter) {
    sigGenes <- c(sigGenes, i)
    outTab <- rbind(
      outTab,
      cbind(
        id = i,
        HR = coxSummary$conf.int[, "exp(coef)"],
        HR.95L = coxSummary$conf.int[, "lower .95"],
        HR.95H = coxSummary$conf.int[, "upper .95"],
        pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
      )
    )
  }
}

outTab <- outTab[order(outTab$HR), ]
write.table(outTab, file = "", sep = "\t", row.names = FALSE, quote = FALSE)
# write.csv(outTab, 'outTab.csv')

uniSigExp <- TME[, sigGenes]
lasso_data <- uniSigExp

## Univariate Cox forest plot 

coxPfilter <- 0.05
outTab <- data.frame()
sigGenes <- c("OS", "OS_event")

for (i in colnames(TME[, 1:37])) {
  cox <- coxph(Surv(OS, OS_event) ~ TME[, i], data = TME)
  coxSummary <- summary(cox)
  coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
  
  sigGenes <- c(sigGenes, i)
  outTab <- rbind(
    outTab,
    cbind(
      id = i,
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
    )
  )
}

outTab <- outTab[order(outTab$HR), ]
outTab$HR <- as.numeric(outTab$HR)
outTab$HR.95L <- as.numeric(outTab$HR.95L)
outTab$HR.95H <- as.numeric(outTab$HR.95H)
outTab$pvalue <- as.numeric(outTab$pvalue)
outTab$pvalue <- ifelse(outTab$pvalue <= 0.05, "P<=0.05", "P>0.05")
outTab$id <- factor(outTab$id, levels = outTab$id)

ggplot(outTab) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_linerange(aes(x = id, ymin = HR.95L, ymax = HR.95H, color = pvalue)) +
  geom_point(aes(x = id, y = HR, color = pvalue)) +
  scale_color_manual(values = c("P<=0.05" = "#d55e00", "P>0.05" = "#0072b2")) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  coord_flip() +
  xlab("Cells in TME") +
  ylab("Hazard Ratio (95% CI)") +
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(file = 'TME_Univariate_Cox_forest_plot .pdf', width = 6, height = 5)


## LASSO Cox regression

set.seed(2000)
x <- as.matrix(lasso_data[, c(3:11)])
y <- as.matrix(Surv(lasso_data$OS, lasso_data$OS_event))

alpha1_fit <- glmnet(x, y, alpha = 1, family = "cox", nlambda = 100)
alpha1.fit.cv <- cv.glmnet(
  x, y,
  type.measure = "deviance",
  alpha = 1,
  family = "cox",
  nfolds = 10
)

print(alpha1.fit.cv)
coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.min)

pdf('Cross_validation_plot.pdf', width = 5, height = 4)
plot(alpha1.fit.cv)
dev.off()

pdf('Coefficient_path_plot.pdf', width = 8, height = 5)
plot_glmnet(alpha1_fit, col = 1:20)
dev.off()

feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.min)))
colnames(feature_all) <- "coff"
feature_opt <- feature_all %>% filter(abs(coff) > 0)
rownames(feature_opt)

feature_ranking <- feature_opt[order(abs(feature_opt$coff), decreasing = TRUE), , drop = FALSE]
print(feature_ranking)

TME_1 <- TME[, rownames(feature_opt)]
TME_2 <- TME[, 38:39]  
lasso <- cbind(TME_1, TME_2)


## Multivariate Cox regression

res.cox <- coxph(
  Surv(OS, OS_event) ~ 
    `CD103+ Resident CD8T` +
    `PDL1+ Mac` +
    `Ki67+ CD4T` +
    `CD11b+ Mac` +
    `CD163hi M2-like RTM` +
    `Undefined CD8T`,
  data = lasso
)

multicox <- summary(res.cox)
print(multicox)


a <- as.data.frame(multicox[["coefficients"]])
b <- as.data.frame(multicox[["conf.int"]])
# write_csv(a, 'a.csv')
# write_csv(b, 'b.csv')


## Multivariate Cox forest plot

forestdata1 <- as.data.frame(multicox[["coefficients"]])
forestdata2 <- as.data.frame(multicox[["conf.int"]])

forestdata1 <- rownames_to_column(forestdata1, var = 'cell')
forestdata2 <- rownames_to_column(forestdata2, var = 'cell')

forestdata1 <- forestdata1[, -c(2, 4, 5)]
forestdata2 <- forestdata2[, c(4, 5)]

forestdata <- cbind(forestdata1, forestdata2)
forestdata <- forestdata[, c(1, 2, 4, 5, 3)]
colnames(forestdata)[3] <- 'lower_95'
colnames(forestdata)[4] <- 'upper_95'
colnames(forestdata)[2] <- 'exp_coef'

forestdata$`Pr(>|z|)` <- ifelse(forestdata$`Pr(>|z|)` <= 0.05, "P<=0.05", "P>0.05")

ggplot(forestdata) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_linerange(aes(x = cell, ymin = lower_95, ymax = upper_95, color = `Pr(>|z|)`)) +
  geom_point(aes(x = cell, y = exp_coef, color = `Pr(>|z|)`)) +
  scale_color_manual(values = c("P<=0.05" = "#d55e00", "P>0.05" = "#0072b2")) +
  scale_y_continuous(limits = c(0, 3.5), expand = c(0, 0)) +
  coord_flip() +
  xlab("Cell Type") +
  ylab("Hazard Ratio (95% CI)") +
  theme_bw()

ggsave(file = 'Multivariate_Cox_forest_plot.pdf', width = 4, height = 2)


## Compute risk score

TME$TME_risk_score <- TME$`PDL1+ Mac` * 0.20489 +
  TME$`CD11b+ Mac` * 0.13681

# Convert OS from days to months
TME$OS <- TME$OS / 30
# write.csv(TME,'155_TME_risk_score.csv')

## Survival analysis 

TME <- read.csv('155_TME_risk_score.csv', row.names = 1, check.names = FALSE)
TME$OS <- TME$OS / 30

survival_time <- TME$OS
survival_event <- TME$OS_event

selected_data <- TME %>%
  select(TME_risk_score) %>%
  mutate(OS = survival_time, OS_event = survival_event) %>%
  na.omit()

best_threshold_surv <- surv_cutpoint(
  selected_data,
  time = "OS",
  event = "OS_event",
  variables = "TME_risk_score",
  minprop = 0.3,
  progressbar = TRUE
)

selected_data <- selected_data %>%
  mutate(
    group = if_else(
      TME_risk_score > best_threshold_surv$cutpoint$cutpoint,
      "high_risk", "low_risk"
    )
  )

selected_data <- selected_data %>%
  mutate(group = factor(group, levels = c('low_risk', 'high_risk'))) %>%
  arrange(group)

# Cox model: group comparison
cox_model_group <- coxph(Surv(OS, OS_event) ~ group, data = selected_data)
summary_cox <- summary(cox_model_group)

hr <- round(summary_cox$coef[1, "exp(coef)"], 2)
ci_lower <- round(summary_cox$conf.int[1, "lower .95"], 2)
ci_upper <- round(summary_cox$conf.int[1, "upper .95"], 2)
pvalue_display <- ifelse(summary_cox$coef[1, "Pr(>|z|)"] < 0.001, "P < 0.001",
                         paste0("P = ", signif(summary_cox$coef[1, "Pr(>|z|)"], 2)))

hr_display <- paste0("HR = ", hr, " (", ci_lower, " - ", ci_upper, ")")
selected_data <- selected_data %>%
  mutate(group = factor(group, levels = c('high_risk', 'low_risk'))) %>%
  arrange(group)

# Kaplan-Meier survival plot
km_fit <- survfit(Surv(OS, OS_event) ~ group, data = selected_data)

plot_object <- ggsurvplot(
  km_fit,
  data = selected_data,
  risk.table = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  pval = paste0(hr_display, "\n", pvalue_display),
  surv.median.line = "hv",
  xlab = 'Follow up times (months)',
  legend.labs = c("high risk", "low risk"),
  risk.table.height = 0.2,
  risk.table.y.text = FALSE,
  ggtheme = theme_survminer() +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold")
    )
)

plot_object

ggsave(
  filename = "155_TME_risk_score_sur_COX.pdf",
  plot = ggarrange(
    plot_object$plot, plot_object$table,
    ncol = 1, nrow = 2,
    heights = c(3, 1)
  ),
  device = "pdf",
  width = 5, height = 6
)

write.csv(TME, '155_TME_risk_score.csv')


## ROC curve evaluation

roc <- roc(TME$OS_event, TME$TME_risk_score)
auc_ci <- ci.auc(roc)
print(auc_ci)

pdf("TME_roc_curve.pdf", width = 6, height = 6)
plot(
  roc,
  col = "#00bfc4",
  print.auc = TRUE,
  auc.polygon = TRUE,
  grid = c(0.1, 0.2),
  grid.col = c("green", "red"),
  max.auc.polygon = TRUE,
  auc.polygon.col = "lightblue",
  print.thres = TRUE
)
par(cex.axis = 4)
par(cex.lab = 4)
par(cex.main = 4)
dev.off()


## Bar plot of risk score coefficients

TME$TME_risk_score <- TME$`PDL1+ Mac` * 0.20489 +
  TME$`CD11b+ Mac` * 0.13681
TME$OS <- TME$OS / 30

df <- data.frame(
  TME = c('PDL1+ Mac', 'CD11b+ Mac'),
  coef = c(0.20489, 0.13681)
)

ggplot(df, aes(x = reorder(TME, coef), y = coef, fill = TME)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CD11b+ Mac" = "#c7e9c0", "PDL1+ Mac" = "#31a354")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.4f", coef)), vjust = -0.5, size = 3) +
  labs(
    x = "Cell Type",
    y = "Coefficient",
    fill = "TME"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

ggsave(filename = "Coefficient_barplot.pdf", width = 4, height = 4)

## ===========================================
###Proteomics riskscore related analysis###
## ===========================================
setwd('D:/BaiduSyncdisk/iCCA/basedprotein')

library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(glmnet)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsignif)

rt <- read.csv('sur_expr.csv', row.names = 1)

# Convert OS to numeric and from days → months
rt$OS <- as.numeric(rt$OS)
rt$OS <- rt$OS / 30
rt <- rt %>% mutate(across(3:ncol(.), ~ log2(. + 1)))
rt_scaled <- as.data.frame(scale(rt[, 3:ncol(rt)]))
rt <- cbind(rt[, 1:2], rt_scaled)

## Univariate Cox regression

coxPfilter = 0.05
outTab = data.frame()
sigGenes = c("OS", "OS_event")

for (i in colnames(rt[, 3:ncol(rt)])) {
  cox <- coxph(Surv(OS, OS_event) ~ rt[, i], data = rt)
  coxSum <- summary(cox)
  pval <- coxSum$coefficients[, "Pr(>|z|)"]
  
  if (pval < coxPfilter) {
    sigGenes <- c(sigGenes, i)
    outTab <- rbind(outTab,
                    cbind(
                      id = i,
                      HR = coxSum$conf.int[, "exp(coef)"],
                      HR.95L = coxSum$conf.int[, "lower .95"],
                      HR.95H = coxSum$conf.int[, "upper .95"],
                      pvalue = pval
                    ))
  }
}

outTab <- outTab[order(outTab$HR), c(1,2,5)]
outTab$logpvalue <- -log10(as.numeric(outTab$pvalue))

uniSigExp <- rt[, sigGenes]
write.csv(uniSigExp, 'uniSigGenes_expr.csv')
# write.csv(outTab, 'Univariate_COX_results.csv')


## LASSO Cox regression

exp_sur <- read.csv('uni_expr_155.csv')
lasso_data <- exp_sur

set.seed(530)
x <- as.matrix(lasso_data[, 4:ncol(lasso_data)])
y <- Surv(lasso_data$OS, lasso_data$OS_event)

fit <- glmnet(x, y, alpha = 1, family = "cox")
cvfit <- cv.glmnet(x, y, type.measure = "deviance", 
                   alpha = 1, family = "cox", nfolds = 5)

pdf('Cross_validation_plot.pdf', width = 5, height = 4)
plot(cvfit)
dev.off()

pdf('Coefficient_path_plot.pdf', width = 5, height = 4)
plot_glmnet(fit, col = 1:20)
dev.off()

feature_all <- as.data.frame(as.matrix(coef(cvfit, s = cvfit$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <- feature_all[abs(feature_all$coff) > 0, , drop = FALSE]

feature_ranking <- feature_opt[order(abs(feature_opt$coff), decreasing = TRUE), ]
write.csv(feature_ranking, "Selected_genes_LASSO.csv")

selected_genes <- c("OS", "OS_event", rownames(feature_ranking))
multidata <- lasso_data[, selected_genes]


## Multivariate Cox regression

res.cox <- coxph(
  Surv(OS, OS_event) ~ SLC2A1 + PLEKHA6 + PAN2 + AGR2 + CLIC3 + SRP14,
  data = multidata
)

multicox <- summary(res.cox)

## Forestplot data
forest1 <- rownames_to_column(as.data.frame(multicox$coefficients), "gene")
forest2 <- rownames_to_column(as.data.frame(multicox$conf.int)[, c(4,5)], "gene")
forest <- cbind(forest1[, c(1,2)], forest2, forest1[, 3])

colnames(forest) <- c("gene","exp_coef","lower_95","upper_95","Pr(>|z|)")
forest$p_col <- ifelse(forest$`Pr(>|z|)` <= 0.05, "P<=0.05", "P>0.05")

ggplot(forest) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_linerange(aes(x = gene, ymin = lower_95, ymax = upper_95, color = p_col)) +
  geom_point(aes(x = gene, y = exp_coef, color = p_col)) +
  scale_color_manual(values = c("P<=0.05" = "#d55e00", "P>0.05" = "#0072b2")) +
  coord_flip() +
  ylab("Hazard Ratio (95% CI)") +
  xlab("Gene Symbol") +
  theme_bw()

ggsave("Multivariate_COX_forest.pdf", width=4, height=2.5)


## Compute protein risk score

exp_sur <- read.csv('protein_risk_score_155.csv', row.names = 1)

exp_sur$protein_risk_score <- 
  exp_sur$PLEKHA6 * -0.2265 +
  exp_sur$PAN2    * -0.3409 +
  exp_sur$CLIC3   *  0.2645 +
  exp_sur$AGR2    *  0.3313


## KM survival analysis

selected_data <- exp_sur %>%
  select(protein_risk_score) %>%
  mutate(OS = exp_sur$OS, OS_event = exp_sur$OS_event) %>%
  na.omit()

cox_model <- coxph(Surv(OS, OS_event) ~ protein_risk_score, selected_data)
selected_data$risk_score <- predict(cox_model, type = "risk")

set.seed(1314)
cut <- surv_cutpoint(selected_data, 
                     time="OS", event="OS_event", 
                     variables="risk_score", minprop=0.3)

selected_data$group <- ifelse(selected_data$risk_score > cut$cutpoint$cutpoint,
                              "high_risk","low_risk")

km_fit <- survfit(Surv(OS, OS_event) ~ group, selected_data)
survd <- survdiff(Surv(OS, OS_event) ~ group, selected_data)
p_val <- 1 - pchisq(survd$chisq, length(survd$n) - 1)

p_display <- ifelse(p_val < 0.001, "<0.001", round(p_val, 3))

plot_object <- ggsurvplot(
  km_fit, data = selected_data,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = paste0("Log-rank\np ", p_display),
  xlab = "Follow up (months)",
  legend.labs = c("high risk","low risk"),
  risk.table.height = 0.2,
  ggtheme = theme_bw()
)

ggsave("protein_survival_curve.pdf",
       ggarrange(plot_object$plot, plot_object$table, ncol=1, heights=c(3,1)),
       width=5, height=6)


## ROC analysis

risk_score <- read.csv('protein_risk_score_155.csv', row.names = 1)

roc_obj <- roc(risk_score$OS_event, risk_score$protein_risk_score)

pdf("protein_roc_curve.pdf", width = 6, height = 6)
plot(roc_obj, col="#00bfc4", print.auc=TRUE, auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres = TRUE)
dev.off()

###clinical risk score related analysis###


## Clean environment and load required packages


setwd('D:/BaiduSyncdisk/iCCA/clinic/')

library(survival)
library(survminer)
library(autoReg)
library(dplyr)
library(flextable)
library(officer)
library(readxl)
library(glmnet)
library(ggplot2)

clinical <- read.csv('clinical_155.csv', row.names = 1)
clinical <- clinical[, -c(3,10)]

clinical$Gender <- factor(clinical$Gender)
clinical$Intrahepatic_Metastasis <- factor(clinical$Intrahepatic_Metastasis)
clinical$Regional_Lymph_Node_Metastasis <- factor(clinical$Regional_Lymph_Node_Metastasis)
clinical$Perineural_Invasion <- factor(clinical$Perineural_Invasion)
clinical$Vascular_Invasion <- factor(clinical$Vascular_Invasion)
clinical$Distal_Metastasis <- factor(clinical$Distal_Metastasis)
clinical$HBV_Status <- factor(clinical$HBV_Status)
clinical$Age <- factor(clinical$Age)
clinical$CA199_cont <- as.numeric(clinical$CA199_cont)
clinical$TNM_Stage <- factor(clinical$TNM_Stage)
clinical$differentiation <- factor(clinical$differentiation)


## Univariate Cox regression

coxPfilter = 0.05
outTab = data.frame()
sigGenes = c("OS", "OS_event")

for (i in colnames(clinical[, 1:11])) {
  cox <- coxph(Surv(OS, OS_event) ~ clinical[, i], data = clinical)
  s <- summary(cox)
  p <- s$coefficients[, "Pr(>|z|)"]
  
  if (p < coxPfilter) {
    sigGenes <- c(sigGenes, i)
    outTab <- rbind(outTab,
                    cbind(id = i,
                          HR = s$conf.int[, "exp(coef)"],
                          HR.95L = s$conf.int[, "lower .95"],
                          HR.95H = s$conf.int[, "upper .95"],
                          pvalue = p))
  }
}

outTab <- outTab[order(outTab$HR), ]
uniclincial <- clinical[, sigGenes]


## LASSO Cox regression

clinical <- uniclincial
clinical$Gender <- as.numeric(factor(clinical$Gender))
clinical$Regional_Lymph_Node_Metastasis <- as.numeric(factor(clinical$Regional_Lymph_Node_Metastasis))
clinical$Intrahepatic_Metastasis <- as.numeric(factor(clinical$Intrahepatic_Metastasis))

lasso_data <- clinical
set.seed(520)

x <- as.matrix(lasso_data[, 3:7])
y <- Surv(lasso_data$OS, lasso_data$OS_event)

fit <- glmnet(x, y, alpha = 1, family = "cox", nlambda = 100)
cvfit <- cv.glmnet(x, y, type.measure = "deviance",
                   alpha = 1, family = "cox", nfolds = 5)

coef(cvfit, s = cvfit$lambda.min)

pdf("Cross_validation_plot.pdf", width = 5, height = 4)
plot(cvfit)
dev.off()

library(plotmo)
pdf("Coefficient_path_plot.pdf", width = 5, height = 5)
plot_glmnet(fit, col = 1:20)
dev.off()


## Clinical risk score

clinical$clinical_risk_score =
  0.1719435 * clinical$Intrahepatic_Metastasis +
  0.5664465 * clinical$Regional_Lymph_Node_Metastasis

clinical$OS <- clinical$OS / 30
write.csv(clinical, "155_clinical_risk_score.csv")

## Survival analysis (KM curves)

selected_data <- clinical %>%
  select(clinical_risk_score) %>%
  mutate(OS = clinical$OS, OS_event = clinical$OS_event) %>%
  na.omit()

cox_model <- coxph(Surv(OS, OS_event) ~ clinical_risk_score, selected_data)
selected_data$risk_score <- predict(cox_model, type = "risk")

set.seed(1314)
cut <- surv_cutpoint(selected_data,
                     time = "OS", event = "OS_event",
                     variables = "risk_score", minprop = 0.3)

selected_data$group <- ifelse(
  selected_data$risk_score > cut$cutpoint$cutpoint,
  "high_risk", "low_risk"
)

selected_data$group <- factor(selected_data$group, levels = c("high_risk", "low_risk"))

km_fit <- survfit(Surv(OS, OS_event) ~ group, selected_data)

plot_object <- ggsurvplot(
  km_fit,
  data = selected_data,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = "HR = 3.22 (2.06 -5.03)\nP < 0.001",
  xlab = "Follow up times (months)",
  legend.labs = c("high risk", "low risk"),
  risk.table.height = 0.2
)


## Univariate Cox forest plot

clinical <- clinical[, -c(3,10)]

outTab <- data.frame()
for (i in colnames(clinical[, 1:11])) {
  cox <- coxph(Surv(OS, OS_event) ~ clinical[, i], data = clinical)
  s <- summary(cox)
  outTab <- rbind(outTab,
                  cbind(id = i,
                        HR = s$conf.int[, "exp(coef)"],
                        HR.95L = s$conf.int[, "lower .95"],
                        HR.95H = s$conf.int[, "upper .95"],
                        pvalue = s$coefficients[, "Pr(>|z|)"]))
}

outTab <- outTab[order(outTab$HR), ]

outTab$HR <- as.numeric(outTab$HR)
outTab$HR.95L <- as.numeric(outTab$HR.95L)
outTab$HR.95H <- as.numeric(outTab$HR.95H)
outTab$pvalue <- ifelse(outTab$pvalue <= 0.05, "P<=0.05", "P>0.05")

ggplot(outTab) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_linerange(aes(x = id, ymin = HR.95L, ymax = HR.95H, color = pvalue)) +
  geom_point(aes(x = id, y = HR, color = pvalue)) +
  scale_color_manual(values = c("P<=0.05" = "#d55e00", "P>0.05" = "#0072b2")) +
  coord_flip() +
  xlab("Clinical Pathology") +
  ylab("Hazard Ratio (95% CI)") +
  theme_bw()

ggsave("Clinical_forest_plot.pdf", width = 6, height = 4)


## Barplot of coefficients

df <- data.frame(
  clinical_pathology = c("Intrahepatic_Metastasis",
                         "Regional_Lymph_Node_Metastasis"),
  coef = c(0.1719435, 0.5664465)
)

ggplot(df, aes(x = reorder(clinical_pathology, coef), y = coef, fill = clinical_pathology)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Intrahepatic_Metastasis" = "#c7e9c0",
                               "Regional_Lymph_Node_Metastasis" = "#31a354")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.4f", coef)), vjust = -0.5, size = 3) +
  labs(x = "Clinical Pathology", y = "Coefficient") +
  theme_classic() +
  theme(axis.text.x = element_blank())

ggsave("Clinical_barplot.pdf", width = 6, height = 4)

## ===========================================
###Spatial multi-modal prognostic scoring system###
## ===========================================

setwd('D:/BaiduSyncdisk/iCCA/155')

library(survival)
library(survminer)
library(dplyr)
library(glmnet)
library(timeROC)
library(ggplot2)
library(tidyverse)
library(flextable)
library(officer)

protein <- read.csv('protein_risk_score_155.csv')
clinical <- read.csv('155_clinical_risk_score.csv')
clinical <- clinical[, c(1, 9)]
TME <- read.csv('155_TME_risk_score.csv')
TME <- TME[, c(1, 41)]

colnames(clinical)[1] <- 'X'
all <- merge(protein, clinical, by = 'X')
all <- merge(all, TME, by = 'X')
colnames(all)[1] <- 'roi_id'

## LASSO Cox model

lasso_data <- all
set.seed(1314)

x <- as.matrix(lasso_data[, 4:6])
y <- Surv(lasso_data$OS, lasso_data$OS_event)

fit <- glmnet(x, y, alpha = 1, family = "cox")
cvfit <- cv.glmnet(x, y, type.measure = "deviance",
                   alpha = 1, family = "cox", nfolds = 5)

print(cvfit)
coef(cvfit, s = cvfit$lambda.min)


## Multi-modal risk score

all$all_risk_score =
  all$protein_risk_score * 1.0032614 +
  all$clinical_risk_score * 1.0409962 +
  all$TME_risk_score * 0.8339307

write.csv(all, '155_all_risk_score.csv')


## Survival analysis

all <- read.csv('155_all_risk_score.csv', row.names = 1)

selected_data <- all %>%
  select(all_risk_score) %>%
  mutate(OS = all$OS,
         OS_event = all$OS_event,
         roi = all$roi_id)

selected_data <- column_to_rownames(selected_data, var = 'roi')
selected_data <- na.omit(selected_data)

cox_model <- coxph(Surv(OS, OS_event) ~ all_risk_score, data = selected_data)
selected_data$risk_score <- predict(cox_model, type = "risk")

hr_val <- round(summary_cox$coef[1, "exp(coef)"], 2)
ci_lower <- round(summary_cox$conf.int[1, "lower .95"], 2)
ci_upper <- round(summary_cox$conf.int[1, "upper .95"], 2)
p_val <- summary_cox$coef[1, "Pr(>|z|)"]
hr_display <- paste0("HR = ", hr_val, " (", ci_lower, " - ", ci_upper, ")")
if (p_val < 0.001) {
  pvalue_display <- "P < 0.001"
} else {
  pvalue_display <- paste0("P = ", round(p_val, 3))
}

set.seed(1314)
cut <- surv_cutpoint(selected_data,
                     time = "OS",
                     event = "OS_event",
                     variables = "risk_score",
                     minprop = 0.3)

selected_data$group <- ifelse(selected_data$risk_score > cut$cutpoint$cutpoint,
                              "high_risk", "low_risk")

selected_data$group <- factor(selected_data$group,
                              levels = c("high_risk", "low_risk"))
selected_data <- selected_data[order(selected_data$group), ]
km_fit <- survfit(Surv(OS, OS_event) ~ group, data = selected_data)

plot_object <- ggsurvplot(
  km_fit,
  data = selected_data,
  risk.table = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  pval = paste0(hr_display, "\n", pvalue_display),
  surv.median.line = "hv",
  xlab = "Follow up times (months)",
  legend.labs = c("high risk", "low risk"),
  risk.table.height = 0.2,
  risk.table.y.text = FALSE,
  ggtheme = theme_survminer()
)

ggsave("155_all_risk_score_sur_COX.pdf",
       ggarrange(plot_object$plot, plot_object$table,
                 ncol = 1, heights = c(3, 1)),
       width = 5, height = 6)


## ROC curve 

library(pROC)

roc_obj <- roc(all$OS_event, all$all_risk_score)

pdf("all_risk_score_roc_curve.pdf", width = 6, height = 6)
plot(roc_obj,
     col = "#00bfc4",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "lightblue",
     print.thres = TRUE)
dev.off()

## Barplot of coefficients

df <- data.frame(
  feature = c("protein_risk_score", "clinical_risk_score", "TME_risk_score"),
  coef = c(1.0032614, 1.0409962, 0.8339307)
)

ggplot(df, aes(x = reorder(feature, coef), y = coef, fill = feature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#fc8d62", "#8da0cb", "#66c2a5")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.4f", coef)), vjust = -0.5, size = 3) +
  labs(x = "", y = "Coefficients", fill = "") +
  theme_classic()

ggsave("all_barplot.pdf", width = 4.5, height = 4)


## ROC comparison

df <- read.csv('155_all_risk_score.csv', row.names = 1)

roc_protein <- roc(df$OS_event ~ df$protein_risk_score)
roc_tme <- roc(df$OS_event ~ df$TME_risk_score)
roc_clinical <- roc(df$OS_event ~ df$clinical_risk_score)
roc_all <- roc(df$OS_event ~ df$all_risk_score)

plot(smooth(roc_protein), col="#fc8d62", lwd=2,
     main="ROC Curve for Risk Scores")
plot(smooth(roc_tme), col="#66c2a5", add=TRUE, lwd=2)
plot(smooth(roc_clinical), col="#8da0cb", add=TRUE, lwd=2)
plot(smooth(roc_all), col="#e78ac3", add=TRUE, lwd=2)

legend("bottomleft",
       legend = c(
         paste0("Proteomic Risk Score: AUC = ", round(auc(roc_protein), 3)),
         paste0("TME Risk Score: AUC = ", round(auc(roc_tme), 3)),
         paste0("Clinical Risk Score: AUC = ", round(auc(roc_clinical), 3)),
         paste0("Multi-modal Risk Score: AUC = ", round(auc(roc_all), 3))
       ),
       col = c("#fc8d62", "#66c2a5", "#8da0cb", "#e78ac3"),
       lty = 1, lwd = 2)

