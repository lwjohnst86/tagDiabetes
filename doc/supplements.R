## ----setup, include=FALSE------------------------------------------------
library(ggplot2)
library(patchwork)
devtools::load_all()
set_options()
knitr::opts_chunk$set(results = "asis", echo = FALSE)

# Captions
sfig <- captioner::captioner(prefix = 'Supplemental Figure')
cite_sf <- pryr::partial(sfig, display = 'cite')

tab <- captioner::captioner(prefix = 'Table')
cite_t <- pryr::partial(tab, display = 'cite')
stab <- captioner::captioner(prefix = 'Supplemental Table')
cite_st <- pryr::partial(stab, display = 'cite')

tab_qic <- stab('qic', 'Comparison of GEE model fitness for variable selection using quasi-likelihood information criteria.')
tab_tagfa <- stab('tagfa', 'Concentration (nmol/mL) and relative percent (mol%) values of TGFA in PROMISE participants at the baseline visit (2004-2006).')
tab_gee_unadj <- stab('gee_unadj', 'Raw estimates and confidence interval values for *time*-adjusted GEE models of the association of the TGFA (mol% and nmol/mL) and total clinically-measured TG with insulin sensitivity and beta-cell function outcomes using the 6 year longitudinal data from the PROMISE cohort. Estimates represent a percent difference in the outcome per SD increase in the FA. P-values were adjusted for the BH false discovery rate, with an asterisk (*) denoting a significant (p<0.05) association.')
tab_gee_adj <- stab('gee_adj', 'Raw estimates and confidence interval values for *fully*-adjusted GEE models of the association of the TGFA (mol% and nmol/mL) and total clinically-measured TG with insulin sensitivity and beta-cell function outcomes using the 6 year longitudinal data from the PROMISE cohort. Variables controlled for were follow-up time, WC, baseline age, ethnicity, sex, ALT, physical activity, and total NEFA. Estimates represent a percent difference in the outcome per SD increase in the FA. P-values were adjusted for the BH false discovery rate, with an asterisk (*) denoting a significant (p<0.05) association.')

fig_consort <- sfig('consort', 'CONSORT diagram of PROMISE participants over the 3 visits.')
fig_dag_is <- sfig('dag_is', 'Directed acyclic graphic output from the DAGitty online software for insulin sensitivity.')
fig_dag_bcf <- sfig('dag_bcf', 'Directed acyclic graphic output from the DAGitty online software for beta-cell function.')
fig_gee_adj_nowaist <- sfig('gee_adj_nowaist', 'Fully-adjusted (without waist size) GEE models of the association of the TGFA (mol% and nmol/mL) and total clinically-measured TG with insulin sensitivity and beta-cell function outcomes using the 6 year longitudinal data from the PROMISE cohort. Variables controlled for were follow-up time, baseline age, ethnicity, sex, ALT, physical activity, and total NEFA. X-axis values represent a percent difference in the outcome per SD increase in the FA. P-values were adjusted for the BH false discovery rate, with the largest dot representing a significant (p<0.05) association.')
fig_pls_loadings <- sfig("pls_loadings", "PLS loadings (or weights) for each of the TGFA. A larger loading indicates a higher contribution to the PLS component score.")

## ----table_qic-----------------------------------------------------------
analyze_qic() %>% 
    table_qic(caption = tab_qic)

## ----table_tagfa---------------------------------------------------------
table_tagfa(project_data, caption = tab_tagfa)

## ----table_gee-----------------------------------------------------------
table_gee_main(gee_results$unadj, caption = tab_gee_unadj)
table_gee_main(gee_results$adj, caption = tab_gee_adj)

## ----fig_gee_adj_nowaist, fig.cap=fig_gee_adj_nowaist--------------------
project_data %>% 
    analyze_gee(covars = covariates[which(covariates != 'Waist')]) %>% 
    dplyr::mutate(p.value = ifelse(p.value <= 0.05, 0.04, p.value)) %>% 
    plot_gee_main() +
    ggplot2::theme(legend.position = "none")

## ----fig_pls_loadings, fig.cap=fig_pls_loadings--------------------------
p1 <- plot_pls_loadings(pls_results$homa2) +
    ggtitle("PLS model loadings for HOMA2-%S")
p2 <- plot_pls_loadings(pls_results$isi) +
    ggtitle("PLS model loadings for ISI")
p1 + p2 + plot_layout(ncol = 1)

## ---- eval=FALSE---------------------------------------------------------
## pls_results$homa2 %>%
##     plot_pls_scores() +
##     theme_classic() +
##     theme(strip.background = element_blank())

