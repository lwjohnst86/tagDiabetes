## ----setup, message=FALSE, echo=FALSE------------------------------------
devtools::load_all('.')
library(ggplot2)
library(patchwork)
library(dplyr)
knitr::opts_chunk$set(echo = FALSE)
ds0yr <- project_data %>% 
    filter(VN == 0)

## Include captions below using `captioner` package
fig <- captioner::captioner(prefix = 'Figure')
cite_f <- pryr::partial(fig, display = 'cite')
sfig <- captioner::captioner(prefix = 'Supplemental Figure S')
cite_sf <- pryr::partial(sfig, display = 'cite')

tab <- captioner::captioner(prefix = 'Table')
cite_t <- pryr::partial(tab, display = 'cite')
stab <- captioner::captioner(prefix = 'Supplemental Table S')
cite_st <- pryr::partial(stab, display = 'cite')
# usage: cite_st('basicChar')

tab_basic <- tab('basic', 'Basic characteristics of PROMISE participants at each of the 3 clinic visits.')

fig_tagfa <- fig('tagfa', 'Distribution of the composition of TGFA in the baseline visit of PROMISE participants (2004-2006). Boxplots represent the median and interquartile range of the FA values.')
fig_heatmap <- fig('heatmap', 'Pearson correlation heatmap of TGFA (nmol/mL) with continuous basic and metabolic characteristics of PROMISE participants from the baseline visit (2004-2006). Darker orange represents a positive correlation; darker blue represents a negative correlation.')
fig_heatmap_mol <- fig('heatmap_mol', 'Pearson correlation heatmap of TGFA (mol%) with continuous basic and metabolic characteristics of PROMISE participants from the baseline visit (2004-2006). Darker orange represents a positive correlation; darker blue represents a negative correlation.')
fig_heatmap_tagfa <- fig('heatmap_tagfa', 'Pearson correlation heatmap of TGFA in the PROMISE participants from the baseline visit (2004-2006). The correlations of FA grouped using heirarchical cluster analysis; FA along the x and y axis are ordered according to this analysis. Darker orange represents a positive correlation; darker blue represents a negative correlation.')
fig_gee_unadj <- fig('gee_unadj', 'Time-adjusted GEE models of the association of the TGFA (mol% and nmol/mL) and total clinically-measured TG with insulin sensitivity and beta-cell function outcomes using the 6 year longitudinal data from the PROMISE cohort. X-axis values represent a percent difference in the outcome per SD increase in the FA. P-values were adjusted for the BH false discovery rate, with the largest dot representing a significant (p<0.05) association.')
fig_gee_adj <- fig('gee_adj', 'Fully-adjusted GEE models of the association of the TGFA (mol% and nmol/mL) and total clinically-measured TG with insulin sensitivity and beta-cell function outcomes using the 6 year longitudinal data from the PROMISE cohort. Variables controlled for were follow-up time, WC, baseline age, ethnicity, sex, ALT, physical activity, and total NEFA. X-axis values represent a percent difference in the outcome per SD increase in the FA. P-values were adjusted for the BH false discovery rate, with the largest dot representing a significant (p<0.05) association.')
fig_pls <- fig('pls', 'Partial least squares (PLS) models showing the clustering of TGFA on insulin sensitivity and beta-cell function measures. The percent explained variance of each component is shown in brackets on each axis. The solid line represents an explained variance of 100% while the dashed line represents an explained variance of 50%. FA between these lines represent variables that strongly explain the underlying structure of the data.')

## ----analyze_others------------------------------------------------------
# Correlation between predicted and observed
pls_homa2_t <- analyze_pls(y = 'lHOMA2_S')
pls_isi_t <- analyze_pls(y = 'lISI')
pls_igi_t <- analyze_pls(y = 'lIGIIR')
pls_issi2_t <- analyze_pls(y = 'lISSI2')

pls_gee <- project_data %>%
    analyze_gee(
        x = c('Comp1', 'Comp2'),
        data_prep = 'pls',
        pls_results = pls_results$isi
        ) %>% 
    dplyr::select(Yterms, Xterms, estimate, p.value) %>% 
    dplyr::filter(Xterms == 'Comp1') %>% 
    dplyr::mutate(p.value = aide::format_pval(p.value),
                  estimate = aide::format_round(estimate))

## ----inline--------------------------------------------------------------
tagfa <- calc_tagfa_percent()
outcome_chg <- calc_outcome_changes()
female <- calc_discr_npct(ds0yr$Sex, 'Female')
euro <- calc_discr_npct(ds0yr$Ethnicity, 'European')
fhd <- calc_discr_npct(ds0yr$FamHistDiab, 'Yes')
gmest0 <- calc_gee_estci(gee_results$unadj)
gmest1 <- calc_gee_estci(gee_results$adj)
fa_cor <- calc_cor(cor_results$conc)
fa_cor_mol <- calc_cor(cor_results$pct)
cor_pls_homa2 <- calc_pred_corr(pls_homa2_t, pls_homa2_t$test_data)
cor_pls_isi <- calc_pred_corr(pls_isi_t, pls_isi_t$test_data)
cor_pls_igi <- calc_pred_corr(pls_igi_t, pls_igi_t$test_data)
cor_pls_issi2 <- calc_pred_corr(pls_issi2_t, pls_issi2_t$test_data)

## ----table_basic---------------------------------------------------------
project_data %>% 
    table_basic(caption = tab_basic)

## ----fig_tagfa, fig.cap=fig_tagfa----------------------------------------
project_data %>% 
    plot_tagfa()

## ----fig_heatmap, fig.cap=fig_heatmap------------------------------------
plot_heatmap(cor_results$conc)

## ----fig_heatmap_mol, fig.cap=fig_heatmap_mol----------------------------
plot_heatmap(cor_results$pct, unit = 'mol%')

## ----fig_heatmap_tagfa, fig.cap=fig_heatmap_tagfa------------------------
plot_heatmap(cor_results$tagfa, text = FALSE)

## ----fig_gee_unadj, fig.cap=fig_gee_unadj--------------------------------
gee_results$unadj %>% 
    mutate(p.value = ifelse(p.value <= 0.05, 0.04, p.value)) %>% 
    plot_gee_main() +
    ggplot2::theme(legend.position = "none")

## ----fig_gee_adj, fig.cap=fig_gee_adj------------------------------------
gee_results$adj %>% 
    mutate(p.value = ifelse(p.value <= 0.05, 0.04, p.value)) %>% 
    plot_gee_main() +
    ggplot2::theme(legend.position = "none")

## ----fig_pls, fig.cap=fig_pls, fig.height=8, fig.width=4-----------------
p1 <- plot_pls_corrcomps(pls_results$homa2) + ggplot2::ggtitle('A: HOMA2-%S')
p2 <- plot_pls_corrcomps(pls_results$isi) + ggplot2::ggtitle('B: ISI')
p1 + p2 + plot_layout(ncol = 1)

