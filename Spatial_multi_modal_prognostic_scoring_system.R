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

###scRNA-seq associated analysis###
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(copykat)

  # plotting
  library(scplotter)
  library(cowplot)
  library(ggpubr)
  library(pheatmap)
  library(colorspace)

  # enrichment
  library(clusterProfiler)
  library(org.Hs.eg.db)

  # GSVA
  library(GSEABase)
  library(GSVA)
  library(BiocParallel)
})

options(stringsAsFactors = FALSE)
set.seed(10086)

############################
# Paths and outputs
############################
root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Input
sce_rds <- file.path(root_dir, "sce.rds")

# Optional inputs
hallmark_gmt <- file.path(root_dir, "ssGSEA", "hallmark.gmt")  # change if needed

# Outputs
out_dir <- file.path(root_dir, "outputs_tumor_epi")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_dir <- file.path(out_dir, "rds")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

############################
# Parameters
############################
marker_genes <- c("PAN2", "PLEKHA6", "CLIC3", "AGR2")
marker_coef  <- c(-0.3409, -0.2265, 0.2645, 0.3313)
names(marker_coef) <- marker_genes

# Epithelial copyKat visualization: Harmony + UMAP
epi_harmony_dims <- 1:20

# Tumor aneuploid subset integration (Seurat v5 IntegrateLayers CCA)
tumor_hvg_nfeatures <- 5000
tumor_integrated_dims <- 1:20

# Tumor clustering
tumor_resolution <- 0.3
tumor_cluster_col <- "RNA_snn_res.0.3"

# Marker calling
markers_min_pct <- 0.3
markers_logfc   <- 1

# Enrichment
go_p_cut <- 0.05
go_q_cut <- 0.05
kegg_p_cut <- 0.05

# ssGSEA
ssgsea_hvg_nfeatures <- 8000

############################
# functions
############################
save_pdf <- function(p, filename, w = 7, h = 6) {
  ggsave(filename = file.path(fig_dir, filename), plot = p, width = w, height = h)
}

safe_read_gmt <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    stop("GMT file not found: ", gmt_path)
  }
  gmt <- clusterProfiler::read.gmt(gmt_path)
  split(gmt$gene, gmt$term)
}

compute_riskscore <- function(seu, genes, coef_vec) {
  exp <- FetchData(seu, vars = genes)
  for (g in genes) {
    if (!g %in% colnames(exp)) stop("Gene not found in object: ", g)
    exp[, g] <- exp[, g] * coef_vec[g]
  }
  rowSums(exp[, genes, drop = FALSE])
}

run_harmony_umap <- function(seu, batch = "orig.ident", dims = 1:20) {
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- RunHarmony(seu, batch)
  seu <- RunUMAP(seu, dims = dims, reduction = "harmony")
  seu
}

run_cca_integrate <- function(seu, split_by = "orig.ident", nfeatures = 5000, dims = 1:20) {
  # Seurat v5 layers workflow
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu[[split_by]][, 1])

  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)

  seu <- IntegrateLayers(
    object = seu,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )

  # re-join layers after integration
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

  seu <- FindNeighbors(seu, dims = dims, reduction = "integrated.cca")
  seu <- RunUMAP(seu, dims = dims, reduction = "integrated.cca")
  seu
}

############################
# Load Seurat object + extract epithelial
############################
if (!file.exists(sce_rds)) stop("sce.rds not found: ", sce_rds)
sce <- readRDS(sce_rds)

if (!"celltype" %in% colnames(sce@meta.data)) {
  stop("Meta column 'celltype' not found in sce.rds.")
}

sce_epi <- subset(sce, subset = celltype == "Epithelial cell")
rm(sce); gc()

# Ensure copyKat exists
if (!"copyKat" %in% colnames(sce_epi@meta.data)) {
  sce_epi$copyKat <- NA
}
sce_epi$copyKat[is.na(sce_epi$copyKat)] <- "diploid"

saveRDS(sce_epi, file.path(rds_dir, "sce_epi.rds"))

############################
# CopyKAT visualization (Harmony UMAP)
############################
sce_epi_h <- run_harmony_umap(sce_epi, batch = "orig.ident", dims = epi_harmony_dims)
Idents(sce_epi_h) <- sce_epi_h$copyKat

p_copykat <- CellDimPlot(sce_epi_h, group_by = "copyKat", reduction = "umap")
save_pdf(p_copykat, "copyKat_epithelial_umap.pdf", w = 7, h = 6)

saveRDS(sce_epi_h, file.path(rds_dir, "sce_epi_harmony_umap.rds"))

############################
# Tumor epithelial (aneuploid) subset + integration
############################
sce_tumor_rds <- file.path(rds_dir, "sce_tumor.rds")

if (file.exists(sce_tumor_rds)) {
  sce_tumor <- readRDS(sce_tumor_rds)
} else {
  sce_tumor <- subset(sce_epi_h, subset = copyKat == "aneuploid")
  saveRDS(sce_tumor, sce_tumor_rds)
}

# Integration via CCA (recommended for multi-dataset tumor epithelial subset)
sce_tumor_int <- run_cca_integrate(
  sce_tumor,
  split_by = "orig.ident",
  nfeatures = tumor_hvg_nfeatures,
  dims = tumor_integrated_dims
)

# Clustering on integrated.cca
sce_tumor_int <- FindClusters(sce_tumor_int, resolution = tumor_resolution)
saveRDS(sce_tumor_int, file.path(rds_dir, "sce_tumor_integrated_cca.rds"))

############################
# Harmony-based alternative 
############################
sce_tumor_h <- sce_tumor
sce_tumor_h <- NormalizeData(sce_tumor_h, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_tumor_h <- FindVariableFeatures(sce_tumor_h)
sce_tumor_h <- ScaleData(sce_tumor_h)
sce_tumor_h <- RunHarmony(sce_tumor_h, "orig.ident")
sce_tumor_h <- FindNeighbors(sce_tumor_h, reduction = "harmony")
sce_tumor_h <- FindClusters(sce_tumor_h, resolution = tumor_resolution)
sce_tumor_h <- RunUMAP(sce_tumor_h, dims = 1:20, reduction = "harmony")
saveRDS(sce_tumor_h, file.path(rds_dir, "sce_tumor_harmony.rds"))

# Use Harmony object as "final tumor" for downstream (matches your original usage)
sce_tumor_final <- sce_tumor_h

############################
#  Rename clusters to C1–C8 
############################
if (!tumor_cluster_col %in% colnames(sce_tumor_final@meta.data)) {
  stop("Expected cluster column not found: ", tumor_cluster_col)
}

Idents(sce_tumor_final) <- sce_tumor_final[[tumor_cluster_col]][, 1]

new.cluster.ids <- c(
  "0" = "C1",
  "1" = "C2",
  "2" = "C3",
  "3" = "C4",
  "4" = "C5",
  "5" = "C6",
  "6" = "C7",
  "7" = "C8"
)

sce_tumor_final <- RenameIdents(sce_tumor_final, new.cluster.ids)
sce_tumor_final[[tumor_cluster_col]] <- as.character(Idents(sce_tumor_final))
saveRDS(sce_tumor_final, file.path(rds_dir, "sce_tumor_C1C8.rds"))

############################
# Marker discovery + heatmap
############################
scRNA.markers <- FindAllMarkers(
  sce_tumor_final,
  only.pos = TRUE,
  min.pct = markers_min_pct,
  logfc.threshold = markers_logfc
)

top_50 <- scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

write.csv(scRNA.markers, file.path(tbl_dir, "tumor_C1C8_all_markers.csv"), row.names = FALSE)
write.csv(top_50, file.path(tbl_dir, "tumor_C1C8_top50_markers.csv"), row.names = FALSE)

# Heatmap top5 (per cluster)
top_5 <- scRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Ensure scale.genes defined before use
scale.genes <- VariableFeatures(sce_tumor_final)
sce_tumor_final <- ScaleData(
  sce_tumor_final,
  features = unique(c(scale.genes, "FAM183A", "MAP3K19", "C1orf194", "PVRL1", "RP1-27K12.2"))
)

p_heat <- DoHeatmap(
  sce_tumor_final,
  features = unique(as.character(top_5$gene)),
  group.by = tumor_cluster_col,
  assay = "RNA"
) + scale_fill_gradientn(colors = c("#ece7f2","#a6bddb","#2b8cbe"))

save_pdf(p_heat, "tumor_C1C8_top_markers_heatmap.pdf", w = 10, h = 5)

############################
#  UMAP + dotplot of PAN2/PLEKHA6/CLIC3/AGR2
############################
p_umap <- CellDimPlot(
  sce_tumor_final,
  group_by = tumor_cluster_col,
  reduction = "umap",
  label = TRUE,
  label_fg = "orange",
  label_bg = "white",
  label_size = 5
)
save_pdf(p_umap, "umap_C1_C8.pdf", w = 7, h = 5)

p_dot <- DotPlot(sce_tumor_final, features = marker_genes) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL)
save_pdf(p_dot, "C1_C8_marker_dotplot.pdf", w = 4, h = 4)

############################
# Riskscore (violin + UMAP overlay)
############################
rs <- compute_riskscore(sce_tumor_final, marker_genes, marker_coef)
sce_tumor_final$riskscore <- rs
saveRDS(sce_tumor_final, file.path(rds_dir, "sce_tumor_C1C8_with_riskscore.rds"))

metadata_subset <- sce_tumor_final@meta.data[, c(tumor_cluster_col, "riskscore")]
col_palette <- c(
  "#a6cee3","#1f78b4","#b2df8a","#33a02c",
  "#fdbf6f","#ff7f00","#fb9a99","#e31a1c","#cab2d6"
)

c2_mean <- mean(metadata_subset$riskscore[metadata_subset[[tumor_cluster_col]] == "C2"], na.rm = TRUE)

p_violin <- ggviolin(
  metadata_subset,
  x = tumor_cluster_col,
  y = "riskscore",
  color = tumor_cluster_col,
  trim = TRUE,
  size = 0.2,
  palette = col_palette,
  scale = "width",
  add = c("mean_sd")
) +
  labs(title = "", x = "", y = "riskscore") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  geom_hline(yintercept = c2_mean, colour = "gray", linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5))

save_pdf(p_violin, "riskscore_C1_C8_violin.pdf", w = 5, h = 4)

# UMAP overlay
umap_df <- as.data.frame(sce_tumor_final@reductions$umap@cell.embeddings)
umap_df$cell <- rownames(umap_df)
rs_df <- sce_tumor_final@meta.data[, c("riskscore"), drop = FALSE] %>%
  rownames_to_column("cell")
plot_df <- left_join(umap_df, rs_df, by = "cell")

p_rs_umap <- ggplot(plot_df) +
  geom_point(aes(x = umap_1, y = umap_2, color = riskscore), size = 0.5, shape = 16) +
  scale_color_viridis_c(option = "inferno", name = "Riskscore") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

save_pdf(p_rs_umap, "Riskscore_umap_tumor_C1_C8.pdf", w = 5, h = 4)

############################
# GO/KEGG enrichment for C2 markers
############################
markers_loose <- FindAllMarkers(
  sce_tumor_final,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.5
)

markers_C2 <- markers_loose %>% filter(cluster == "C2")

ids <- bitr(markers_C2$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)
markers_C2_ids <- merge(markers_C2, ids, by.x = "gene", by.y = "SYMBOL")
gene_list <- unique(markers_C2_ids$ENTREZID)

go_res <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = go_p_cut,
  qvalueCutoff = go_q_cut
)
go_tbl <- as.data.frame(go_res)
write.csv(go_tbl, file.path(tbl_dir, "GO_C2_full.csv"), row.names = FALSE)

kegg_res <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",
  pvalueCutoff = kegg_p_cut
)
kegg_tbl <- as.data.frame(kegg_res)
write.csv(kegg_tbl, file.path(tbl_dir, "KEGG_C2_full.csv"), row.names = FALSE)

############################
# Hallmark ssGSEA (GSVA) + heatmap + correlation with riskscore (C2)
############################
if (file.exists(hallmark_gmt)) {
  genesets_hall <- safe_read_gmt(hallmark_gmt)

  sce_tumor_final <- FindVariableFeatures(sce_tumor_final, nfeatures = ssgsea_hvg_nfeatures)
  scale.genes2 <- VariableFeatures(sce_tumor_final)

  expr <- as.matrix(GetAssayData(sce_tumor_final, layer = "data", assay = "RNA")[scale.genes2, ])

  params <- gsvaParam(
    expr = expr,
    genesets_hall,
    minSize = 3,
    maxSize = Inf,
    kcdf = "Gaussian",
    tau = 1,
    maxDiff = TRUE,
    absRanking = FALSE
  )

  gsva_result <- gsva(
    params,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
  )

  ssgsea_scores <- as.data.frame(t(gsva_result))
  sce_tumor_final <- AddMetaData(sce_tumor_final, ssgsea_scores)

  cluster_info_df <- data.frame(
    cell_name = colnames(sce_tumor_final),
    cluster = sce_tumor_final@meta.data[[tumor_cluster_col]],
    stringsAsFactors = FALSE
  )

  ssgsea_scores$cell_name <- rownames(ssgsea_scores)
  ssgsea_scores <- merge(ssgsea_scores, cluster_info_df, by = "cell_name")
  rownames(ssgsea_scores) <- ssgsea_scores$cell_name
  ssgsea_scores$cell_name <- NULL

  avg_ssgsea_scores <- ssgsea_scores %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames("cluster")

  pdf(file.path(fig_dir, "hallmark_C1_C8_heatmap.pdf"), width = 8, height = 12)
  pheatmap(
    t(avg_ssgsea_scores),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = "white",
    color = colorRampPalette(c("#2ca02c", "white", "#ff7f0e"))(100)
  )
  dev.off()

  # Correlation with riskscore within C2
  ssgsea_C2 <- ssgsea_scores %>% filter(cluster == "C2")
  # drop cluster column
  ssgsea_C2_mat <- ssgsea_C2 %>% select(-cluster)

  # align riskscore
  rs_C2 <- sce_tumor_final@meta.data[sce_tumor_final@meta.data[[tumor_cluster_col]] == "C2", "riskscore"]
  n_cor <- cbind(riskscore = rs_C2, ssgsea_C2_mat)

  cor_matrix <- cor(n_cor %>% select(-riskscore), n_cor$riskscore, method = "pearson")

  pdf(file.path(fig_dir, "Pearson_Correlation_C2_heatmap.pdf"), width = 6, height = 12)
  pheatmap(
    as.matrix(cor_matrix),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    border_color = "white",
    number_format = "%.2f",
    color = colorRampPalette(c("#2ca02c", "white", "#ff7f0e"))(50),
    main = "Pearson correlation with riskscore (C2)",
    fontsize_number = 10,
    fontsize = 12
  )
  dev.off()

  # Example scatter panels (edit pathways if needed)
  pick_paths <- c(
    "ESTROGEN_RESPONSE_LATE",
    "ESTROGEN_RESPONSE_EARLY",
    "KRAS_SIGNALING_UP",
    "TNFA_SIGNALING_VIA_NFKB"
  )
  pick_paths <- pick_paths[pick_paths %in% colnames(n_cor)]

  plist <- list()
  colors_sc <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
  for (i in seq_along(pick_paths)) {
    y_var <- pick_paths[i]
    p <- ggscatter(
      n_cor, x = "riskscore", y = y_var,
      add = "reg.line", conf.int = TRUE,
      add.params = list(fill = "lightgray"),
      color = colors_sc[min(i, length(colors_sc))],
      alpha = 0.1,
      size = 0.7
    ) +
      stat_cor(method = "pearson") +
      labs(title = y_var)
    plist[[i]] <- p
  }

  if (length(plist) > 0) {
    p_arr <- ggarrange(plotlist = plist, ncol = 2, nrow = ceiling(length(plist) / 2))
    save_pdf(p_arr, "pearson_selected_pathways_C2.pdf", w = 6, h = 6)
  }

} else {
  message("Hallmark GMT not found, skipping ssGSEA/GSVA: ", hallmark_gmt)
}

############################
#  Coef contribution (C2)
############################
dir.create(file.path(out_dir, "contribution"), showWarnings = FALSE, recursive = TRUE)
contrib_dir <- file.path(out_dir, "contribution")

scRNA_C2 <- subset(sce_tumor_final, subset = .data[[tumor_cluster_col]] == "C2")
exp_C2 <- FetchData(scRNA_C2, vars = marker_genes)

for (g in marker_genes) exp_C2[, g] <- exp_C2[, g] * marker_coef[g]
gene_sum <- colSums(exp_C2)

df_weight <- data.frame(marker = marker_genes, sum = as.numeric(gene_sum))
p_w <- ggplot(df_weight, aes(x = marker, y = sum, fill = marker)) +
  geom_col(position = "fill") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(contrib_dir, "C2_coef_weight.pdf"), plot = p_w, width = 3, height = 3)

df_stack <- df_weight
df_stack$sum_signed <- df_stack$sum
df_stack$sum_signed[df_stack$marker %in% names(marker_coef[marker_coef < 0])] <- -abs(df_stack$sum_signed[df_stack$marker %in% names(marker_coef[marker_coef < 0])])

p_s <- ggplot(df_stack, aes(x = marker, y = sum_signed, fill = marker)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(contrib_dir, "C2_coef_exp.pdf"), plot = p_s, width = 3, height = 3)

