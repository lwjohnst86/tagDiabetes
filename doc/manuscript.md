
<!--

Outline:

- Primary analysis:
    - Heatmap
    - GEE
        - Include # of both uncorrected and corrected (FDR) significant associations
        - QIC modeling
        - DAG development
        - Look over the entire model (see which variables are significant, etc)
    - Interaction (manhattan plot)
    - Interaction plots
    - PLS with ISI and HOMA-IR
        - CV for determining number of comps
        - Loadings (of TAGFA)
        - Prediction (of predicted vs real)
        - Use extracted PLS components in GEE
        - Train and test (CV) to see how it predicts in a new dataset (ie. a future visit)
- Secondary analysis:
    - PLS-DA (dysgly)
        - Train and test (CV)
        - Loadings (of TAGFA) 
        - Groupings
        - Prediction (misclassification)
    - AROC
        - dysgly with clinical TAG vs the 'best' fatty acid to compare to
    
-->

# Title page

Title:

Author: Luke W. Johnston (1), MSc; Stewart B. Harris (2), MD; Ravi Retnakaran
(3,4), MD; Bernard Zinman (3,4), MD; Adria Giacca (5), PhD; Zhen Liu (1), PhD;
Richard P. Bazinet (1), PhD; and Anthony J. Hanley (1,6), PhD

<!-- During time of study. -->
Affiliation:

(1) Department of Nutritional Sciences, University of Toronto, Toronto, Ontario, Canada.
(2) Centre for Studies in Family Medicine, University of Western Ontario, London, Ontario, Canada.
(3) Division of Endocrinology, University of Toronto, Toronto, ON, Canada.
(4) Lunenfeld Tanenbaum Research Institute, Mount Sinai Hospital, Toronto, Ontario, Canada.
(5) Department of Physiology, University of Toronto, Toronto, Ontario, Canada.
(6) Dalla Lana School of Public Health, University of Toronto, Toronto, Ontario, Canada.

Corresponding author:

- Name: Anthony J. Hanley
- Current address:  
    Department of Nutritional Sciences  
    Faculty of Medicine  
    University of Toronto  
    FitzGerald Building, 150 College Street, Room 341  
    Toronto, ON, Canada, M5S 3E2
- Phone number: 416.978.3616
- Fax number:
- Email: anthony.hanley@utoronto.ca

Disclaimers:

Funding support:

# Abstract

# Background

# Subjects and Methods


```r
# Only set if the Rmd file is not in the parent directory (ie. 'projectname/')
devtools::load_all('.')
```

```
## Loading tagDiabetes
```

```r
ls(getNamespace('tagDiabetes'), all.names = TRUE)
```

```
##  [1] "analyze_corr"                       
##  [2] "analyze_gee"                        
##  [3] "analyze_pls"                        
##  [4] "analyze_plsda"                      
##  [5] "calculate_conversion_dysgly"        
##  [6] "calculate_factor_npercent"          
##  [7] "calculate_gee_magnitude"            
##  [8] "calculate_ngroups_lcmm"             
##  [9] "calculate_outcomes_pct_change"      
## [10] "calculate_percent_nefa_contribution"
## [11] "calculate_plsda_misclass"           
## [12] ".__DEVTOOLS__"                      
## [13] "example_plsda_results"              
## [14] "extract_gee_estimateCI"             
## [15] "fetch_data"                         
## [16] "get_dysglycemia_data"               
## [17] "get_gee_data"                       
## [18] "get_pls_data"                       
## [19] "graph_theme"                        
## [20] "load_data"                          
## [21] ".__NAMESPACE__."                    
## [22] ".packageName"                       
## [23] "plot_gee_main"                      
## [24] "plot_heatmap"                       
## [25] "plot_plsda_grouping"                
## [26] "plot_plsda_loadings"                
## [27] "plot_tag_distribution"              
## [28] "renaming_fa"                        
## [29] "renaming_fraction"                  
## [30] "renaming_list"                      
## [31] "renaming_outcomes"                  
## [32] "renaming_table_rows"                
## [33] ".__S3MethodsTable__."               
## [34] "set_options"                        
## [35] "table_basic"                        
## [36] "table_corr"                         
## [37] "table_distribution"                 
## [38] "table_gee_main"                     
## [39] "tidy_gee_results"                   
## [40] "trim_ws"
```

```r
knitr::opts_knit$set(root.dir = '../')
knitr::opts_chunk$set(eval = FALSE)
```


```r
devtools::load_all('.')
run_setup()
ds <- load_data()
```


```r
## Include captions below using `captioner` package
fig <- captioner::captioner(prefix = 'Figure')
cite_f <- pryr::partial(fig, display = 'cite')
sfig <- captioner::captioner(prefix = 'Supplemental Figure')
cite_sf <- pryr::partial(sfig, display = 'cite')
tab <- captioner::captioner(prefix = 'Table')
cite_t <- pryr::partial(tab, display = 'cite')
stab <- captioner::captioner(prefix = 'Supplemental Table')
cite_st <- pryr::partial(stab, display = 'cite')
# usage: cite_st('basicChar')

fig_gee <- fig('gee', 'caption here')
```


```r
outcomes <- c('linvHOMA', 'lISI', 'lIGIIR', 'lISSI2')
tg_conc <- grep('^tg\\d+', names(ds), value = TRUE)
tg_pct <- grep('^pct_tg\\d+', names(ds), value = TRUE)
tg_total <- c('TotalTG', 'TAG') # try with lTAG too.. but exponentiating gives diff interpret..
covariates <- c('VN', 'BaseAge', 'Sex', 'Ethnicity', 'Waist')
```


```r
gee_df <- get_gee_data(ds)
pls_df_homa <- get_pls_data(ds, 'linvHOMA')
pls_df_isi <- get_pls_data(ds, 'lISI')
```


```r
gee_results <-
    dplyr::bind_rows(
        analyze_gee(gee_df, outcomes, tg_conc, covariates, 'nmol/mL'),
        analyze_gee(gee_df, outcomes, tg_pct, covariates, 'mol%'),
        analyze_gee(gee_df, outcomes, tg_total, covariates, 'Totals')
    ) %>%
    tidy_gee_results()

gee_results_int <- 
    dplyr::bind_rows(
        analyze_gee(gee_df, outcomes, tg_conc, covariates, 'nmol/mL', intvar = 'VN'),
        analyze_gee(gee_df, outcomes, tg_pct, covariates, 'mol%', intvar = 'VN'),
        analyze_gee(gee_df, outcomes, tg_total, covariates, 'Totals', intvar = 'VN')
    ) %>%
    tidy_gee_results()
```


```r
pls_results <- list()
pls_results$homa$train <- analyze_pls(pls_df_homa$train, 'linvHOMA')
pls_results$homa$test <- analyze_pls(pls_df_homa$test, 'linvHOMA')
pls_results$isi$train <- analyze_pls(pls_df_isi$train, 'lISI')
pls_results$isi$test <- analyze_pls(pls_df_isi$test, 'lISI')
```


```r
cor_results <- analyze_corr(ds)
```


# Results


## GEE

Main effect

Interaction effect

## PLS

Accuracy of predictive ability of PLS components with real measure {{use test maybe?}}.
Correlation of pred vs real in test set.

## PLS-DA

Misclassification of dysgly {{using test?}}, misclass by group too ('yes' vs 'no')

## AROC

AROC of clinical TAG vs TAGFA palmitic acid etc.

# Discussion

## Acknowledgements

# Tables and figures


```r
table_basic(ds, '')
```


```r
plot_heatmap(cor_results)
```


```r
#table_corr(cor_results)
```


```r
plot_gee_main(gee_results)
```


```r
table_gee_main(gee_results, '')
```


```r
seer::view_interaction(gee_results_int, 
                       groups = 'unit~Yterms',
                       ylab = 'Triacylglycerol fatty acids') +
        graph_theme()
```


```r
# plot_interactions(ds, fa, outcome)
# plot_interactions(ds, fa, outcome)
# plot_interactions(ds, fa, outcome)
```


```r
# Or as a table?
plot_pls_cv_ncomp(pls_results$homa$train)
plot_pls_cv_ncomp(pls_results_homa)
```


```r
# As a single figure
seer::view_pls_xloadings(pls_results$homa$train) +
    graph_theme()
seer::view_pls_xloadings(pls_results$homa$test) +
    graph_theme()
seer::view_pls_xloadings(pls_results$isi$train) +
    graph_theme()
seer::view_pls_xloadings(pls_results$isi$test) +
    graph_theme()
```


```r
# maybe not a figure, maybe just a correlation?
predplot(pls_results$isi$train, ncomp = 1:2)
predplot(pls_results$homa$train, ncomp = 1:2)
plot_pls_predict(pls_results$isi$train)
plot_pls_predict(pls_results_homa)
plot_pls_predict(test_set, pls_predict_set)
```


```r
plot_gee_main(gee_results_plscomps)
```


```r
table_gee_main(gee_results_plscomps)
```


```r
# can't remember if this is possible...
plot_plsda_cv_ncomp(plsda_results)
```


```r
plot_plsda_loadings(plsda_results)
```


```r
plot_plsda_groupings(plsda_results)
```


```r
plot_aroc(aroc_results) #?
```

# Supplemental Methods



# References
