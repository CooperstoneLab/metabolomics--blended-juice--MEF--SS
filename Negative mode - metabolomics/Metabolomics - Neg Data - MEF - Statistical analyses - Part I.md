# Metabolomic - Statistics - Negative mode - MEF - Part I
Giovana Domeneghini Mercali

<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/htmlwidgets-1.6.2/htmlwidgets.js"></script>
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/plotly-binding-4.10.2/plotly.js"></script>
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/setprototypeof-0.1/setprototypeof.js"></script>
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/typedarray-0.1/typedarray.min.js"></script>
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/jquery-3.5.1/jquery.min.js"></script>
<link href="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/crosstalk-1.2.0/css/crosstalk.min.css" rel="stylesheet" />
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/crosstalk-1.2.0/js/crosstalk.min.js"></script>
<link href="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/plotly-htmlwidgets-css-2.11.1/plotly-htmlwidgets.css" rel="stylesheet" />
<script src="Metabolomics - Neg Data - MEF - Statistical analyses - Part I_files/libs/plotly-main-2.11.1/plotly-latest.min.js"></script>


### Load libraries

``` r
library(factoextra) # visualizing PCA results
library(glue) # for easy pasting
library(plotly) # quick interactive plots
library(proxyC) # more efficient large matrices
library(data.table) # easy transposing
library(janitor) # for cleaning names and checking duplicates
library(notame) # for collapsing ions coming from the same metabolite
library(doParallel) # for parallelizing notame specifically
library(patchwork) # for making multi-panel plots
library(readxl)
library(svglite)
library(pheatmap) # for heínsatmaps
library(ggridges) # for ridgeline plots
library(ggdist) # for nice dotplots
library(rstatix) # for pipeable stats testing
library(agricolae) # for posthoc tests 
library(ggpubr) # extension for adding stats to plots


# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

## Read in data

Negative mode.

``` r
# read in metabolomics data
feature_table_long <- read_csv("features_clustered_long__neg_MEF.csv",
                  trim_ws = TRUE,
                  na = "0") # read in 0s to be NAs.
```

    Rows: 75405 Columns: 8
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (4): sample_name, treatment, sample_or_qc, mz_rt
    dbl (4): run_order, mz, rt, rel_abund

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data
metadata <- read_csv("metadata_neg_mef.csv",
                  trim_ws = TRUE)
```

    Rows: 67 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (3): sample_name, short_sample_name, treatment
    dbl (3): run_order, soluble_solids, mass

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data
annotation <- read_csv("Annotation_neg.csv", 
                  locale=locale(encoding = "latin1"))
```

    Rows: 38 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (2): Compound, mz_rt

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Take a quick look at our data.

``` r
# look at first 5 rows, first 5 columns 
feature_table_long[1:5,1:5]
```

    # A tibble: 5 × 5
      sample_name treatment run_order sample_or_qc    mz
      <chr>       <chr>         <dbl> <chr>        <dbl>
    1 MT-2        MT               16 Sample        100.
    2 MT-2        MT               16 Sample        101.
    3 MT-2        MT               16 Sample        101.
    4 MT-2        MT               16 Sample        101.
    5 MT-2        MT               16 Sample        101.

``` r
# look at first 5 rows, all columns 
metadata[1:5,]
```

    # A tibble: 5 × 6
      sample_name        short_sample_name treatment run_order soluble_solids  mass
      <chr>              <chr>             <chr>         <dbl>          <dbl> <dbl>
    1 FreshJuice2_neg_54 FJ-2              FJ               54           10.3 0.428
    2 40TT_10M6_neg_72   MT-6              MT               72           10.2 0.426
    3 40TT4_neg_49       40TT-4            40TT             49            9.9 0.425
    4 40TT1_neg_53       40TT-1            40TT             53           10.1 0.434
    5 40TT6_neg_69       40TT-6            40TT             69           10   0.429

``` r
# check dimensions
dim(feature_table_long)
```

    [1] 75405     8

``` r
dim(metadata)
```

    [1] 67  6

First we are going to convert the type of some of the columns to match
what we want (e.g., run order converted to numeric, treatment to
factor).

``` r
tibble(feature_table_long)
```

    # A tibble: 75,405 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <chr>         <dbl> <chr>        <dbl> <dbl> <chr>          <dbl>
     1 MT-2        MT               16 Sample        100. 0.898 100.0402_…   114435.
     2 MT-2        MT               16 Sample        101. 0.338 100.9439_…    17707.
     3 MT-2        MT               16 Sample        101. 0.321 101.0222_…   776542.
     4 MT-2        MT               16 Sample        101. 0.478 101.0238_…  1110223.
     5 MT-2        MT               16 Sample        101. 0.477 101.031_0…    86064.
     6 MT-2        MT               16 Sample        103. 0.444 103.003_0…    40976.
     7 MT-2        MT               16 Sample        103. 0.362 103.0384_…   135118.
     8 MT-2        MT               16 Sample        103. 0.354 103.0396_…   172927.
     9 MT-2        MT               16 Sample        103. 0.620 103.0397_…   103613.
    10 MT-2        MT               16 Sample        103. 0.856 103.0399_…   309728.
    # ℹ 75,395 more rows

``` r
# make run_order numeric
feature_table_long$run_order <- as.numeric(feature_table_long$run_order)

# make treatment and sample_or_qc a factor (i.e., categorical)
feature_table_long$treatment <- as.factor(feature_table_long$treatment)
feature_table_long$sample_or_qc <- as.factor(feature_table_long$sample_or_qc)

# did it work?
tibble(feature_table_long)
```

    # A tibble: 75,405 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <fct>         <dbl> <fct>        <dbl> <dbl> <chr>          <dbl>
     1 MT-2        MT               16 Sample        100. 0.898 100.0402_…   114435.
     2 MT-2        MT               16 Sample        101. 0.338 100.9439_…    17707.
     3 MT-2        MT               16 Sample        101. 0.321 101.0222_…   776542.
     4 MT-2        MT               16 Sample        101. 0.478 101.0238_…  1110223.
     5 MT-2        MT               16 Sample        101. 0.477 101.031_0…    86064.
     6 MT-2        MT               16 Sample        103. 0.444 103.003_0…    40976.
     7 MT-2        MT               16 Sample        103. 0.362 103.0384_…   135118.
     8 MT-2        MT               16 Sample        103. 0.354 103.0396_…   172927.
     9 MT-2        MT               16 Sample        103. 0.620 103.0397_…   103613.
    10 MT-2        MT               16 Sample        103. 0.856 103.0399_…   309728.
    # ℹ 75,395 more rows

## Data quality

#### Log2 transformed

We will log2 transform our data.

``` r
feature_table_long_log2 <- feature_table_long %>%
  mutate(rel_abund = log2(rel_abund))
```

And then plot.

``` r
boxplot_log2_transformed <- feature_table_long_log2 %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5",
                                "QC" = "#F6A900")) +
  labs(y = "Relative abundance", x = "Sample", 
  title = "LC-MS (-) feature abundances by Sample",
       subtitle = "Data is log2 transformed")
```

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
boxplot_log2_transformed_run_order <-feature_table_long_log2 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5",
                                "QC" = "#F6A900")) +
  labs(y = "Relative abundance", x = "Sample")
```

``` r
boxplot <- boxplot_log2_transformed / boxplot_log2_transformed_run_order + plot_annotation(tag_levels = c("A"))

ggsave(plot = boxplot,
       filename = "boxplot_final_neg.svg")
```

    Saving 7 x 5 in image

Let’s write out that data as our final feature table used for the rest
of the analysis.

``` r
feature_table_wide_log2 <- feature_table_long_log2 %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = mz_rt,
              values_from = rel_abund)

write_csv(feature_table_wide_log2,
          file = "feature_table_wide_log2_neg.csv")
```

## PCAs all data

### With QCs

``` r
pca_qc <- prcomp(feature_table_wide_log2[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_qc)
```

    Importance of components:
                               PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     24.2071 10.6498 7.02727 6.03875 3.36391 2.48708 2.32426
    Proportion of Variance  0.6671  0.1291 0.05622 0.04152 0.01288 0.00704 0.00615
    Cumulative Proportion   0.6671  0.7962 0.85245 0.89396 0.90685 0.91389 0.92004
                               PC8     PC9    PC10    PC11    PC12    PC13    PC14
    Standard deviation     2.24795 2.05186 2.02302 1.93214 1.92363 1.86205 1.82313
    Proportion of Variance 0.00575 0.00479 0.00466 0.00425 0.00421 0.00395 0.00378
    Cumulative Proportion  0.92579 0.93058 0.93524 0.93949 0.94371 0.94765 0.95144
                              PC15   PC16   PC17    PC18    PC19    PC20    PC21
    Standard deviation     1.80001 1.7531 1.7281 1.68970 1.68479 1.64350 1.62994
    Proportion of Variance 0.00369 0.0035 0.0034 0.00325 0.00323 0.00308 0.00302
    Cumulative Proportion  0.95513 0.9586 0.9620 0.96528 0.96851 0.97158 0.97461
                              PC22    PC23    PC24    PC25    PC26    PC27   PC28
    Standard deviation     1.57885 1.52772 1.49466 1.46354 1.41889 1.40623 1.3917
    Proportion of Variance 0.00284 0.00266 0.00254 0.00244 0.00229 0.00225 0.0022
    Cumulative Proportion  0.97744 0.98010 0.98264 0.98508 0.98737 0.98963 0.9918
                              PC29    PC30    PC31   PC32      PC33
    Standard deviation     1.37087 1.35399 1.33812 1.2932 4.654e-14
    Proportion of Variance 0.00214 0.00209 0.00204 0.0019 0.000e+00
    Cumulative Proportion  0.99397 0.99606 0.99810 1.0000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_qc <- summary(pca_qc)$importance %>%
  as.data.frame()

head(importance_qc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     24.20711 10.64978 7.027267 6.038753 3.363915 2.487083
    Proportion of Variance  0.66711  0.12912 0.056220 0.041520 0.012880 0.007040
    Cumulative Proportion   0.66711  0.79623 0.852450 0.893960 0.906850 0.913890
                                PC7      PC8      PC9    PC10     PC11     PC12
    Standard deviation     2.324255 2.247948 2.051856 2.02302 1.932138 1.923632
    Proportion of Variance 0.006150 0.005750 0.004790 0.00466 0.004250 0.004210
    Cumulative Proportion  0.920040 0.925790 0.930580 0.93524 0.939490 0.943710
                               PC13     PC14     PC15     PC16     PC17    PC18
    Standard deviation     1.862046 1.823135 1.800007 1.753076 1.728062 1.68970
    Proportion of Variance 0.003950 0.003780 0.003690 0.003500 0.003400 0.00325
    Cumulative Proportion  0.947650 0.951440 0.955130 0.958630 0.962020 0.96528
                               PC19     PC20     PC21     PC22     PC23     PC24
    Standard deviation     1.684794 1.643503 1.629944 1.578848 1.527719 1.494656
    Proportion of Variance 0.003230 0.003080 0.003020 0.002840 0.002660 0.002540
    Cumulative Proportion  0.968510 0.971580 0.974610 0.977440 0.980100 0.982640
                               PC25     PC26    PC27    PC28     PC29    PC30
    Standard deviation     1.463539 1.418895 1.40623 1.39167 1.370874 1.35399
    Proportion of Variance 0.002440 0.002290 0.00225 0.00220 0.002140 0.00209
    Cumulative Proportion  0.985080 0.987370 0.98963 0.99183 0.993970 0.99606
                               PC31     PC32         PC33
    Standard deviation     1.338117 1.293193 4.654207e-14
    Proportion of Variance 0.002040 0.001900 0.000000e+00
    Cumulative Proportion  0.998100 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-12-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-13-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_qc <- as.data.frame(pca_qc$x)

# bind meta-data
scores_qc <- bind_cols(feature_table_wide_log2[,1:4], # first 4 columns
                       scores_raw_qc)
```

Plot.

``` r
# create objects indicating percent variance explained by PC1 and PC2
PC1_percent_qc <- round((importance_qc[2,1])*100, # index 2nd row, 1st column, times 100
                         1) # round to 1 decimal
PC2_percent_qc <- round((importance_qc[2,2])*100, 1) 

# plot
# aes(text) is for setting tooltip with plotly later to indicate hover text
(scores_qc_plot <- scores_qc %>%
  ggplot(aes(x = PC1, y = PC2, fill = treatment, text = glue("Sample: {sample_name},
                                                           Treatment: {treatment}"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 3, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5",
                                "QC" = "#F6A900")) +
  labs(x = glue("PC1: {PC1_percent_qc}%"), 
       y = glue("PC2: {PC2_percent_qc}%"), 
       title = "", tag = "A"))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-15-1.png)

Then make your scores plot interactive so you can see who is who.

``` r
ggplotly(scores_qc_plot, tooltip = "text")
```

``` r
ggsave(plot = scores_qc_plot,
       filename = "pca_qc_final_neg.jpeg", width = 140, height = 120, units = "mm")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-18-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

### Without QCs

``` r
feature_table_wide_log2_noqc <- feature_table_wide_log2 %>%
  filter(sample_or_qc == "Sample")


pca_noqc <- prcomp(feature_table_wide_log2_noqc[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_noqc)
```

    Importance of components:
                               PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     30.5681 12.2814 8.57791 4.53679 3.26833 2.90375 2.83386
    Proportion of Variance  0.7343  0.1185 0.05782 0.01617 0.00839 0.00663 0.00631
    Cumulative Proportion   0.7343  0.8528 0.91067 0.92685 0.93524 0.94187 0.94818
                               PC8     PC9    PC10    PC11    PC12    PC13  PC14
    Standard deviation     2.61503 2.50201 2.41632 2.40537 2.34358 2.28754 2.256
    Proportion of Variance 0.00537 0.00492 0.00459 0.00455 0.00432 0.00411 0.004
    Cumulative Proportion  0.95355 0.95847 0.96306 0.96761 0.97193 0.97604 0.980
                              PC15    PC16    PC17    PC18    PC19    PC20
    Standard deviation     2.16276 2.12959 2.12624 2.06323 1.98795 1.86030
    Proportion of Variance 0.00368 0.00356 0.00355 0.00335 0.00311 0.00272
    Cumulative Proportion  0.98371 0.98728 0.99083 0.99417 0.99728 1.00000
                                PC21
    Standard deviation     4.574e-14
    Proportion of Variance 0.000e+00
    Cumulative Proportion  1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_noqc <- summary(pca_noqc)$importance %>%
  as.data.frame()

head(importance_noqc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     30.56805 12.28144 8.577913 4.536787 3.268331 2.903749
    Proportion of Variance  0.73432  0.11853 0.057820 0.016170 0.008390 0.006630
    Cumulative Proportion   0.73432  0.85285 0.910670 0.926850 0.935240 0.941870
                                PC7      PC8      PC9     PC10     PC11     PC12
    Standard deviation     2.833862 2.615026 2.502006 2.416319 2.405374 2.343578
    Proportion of Variance 0.006310 0.005370 0.004920 0.004590 0.004550 0.004320
    Cumulative Proportion  0.948180 0.953550 0.958470 0.963060 0.967610 0.971930
                               PC13     PC14     PC15    PC16     PC17     PC18
    Standard deviation     2.287541 2.255655 2.162762 2.12959 2.126238 2.063228
    Proportion of Variance 0.004110 0.004000 0.003680 0.00356 0.003550 0.003350
    Cumulative Proportion  0.976040 0.980040 0.983710 0.98728 0.990830 0.994170
                              PC19     PC20         PC21
    Standard deviation     1.98795 1.860302 4.574379e-14
    Proportion of Variance 0.00311 0.002720 0.000000e+00
    Cumulative Proportion  0.99728 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-21-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-22-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_noqc <- as.data.frame(pca_noqc$x)

# bind meta-data
scores_noqc <- bind_cols(feature_table_wide_log2_noqc[,1:4], # metadata
                         scores_raw_noqc)
```

Plot.

``` r
# create objects indicating percent variance explained by PC1 and PC2
PC1_percent_noqc <- round((importance_noqc[2,1])*100, # index 2nd row, 1st column, times 100
                         1) # round to 1 decimal
PC2_percent_noqc <- round((importance_noqc[2,2])*100, 1) 

# plot
# aes(text) is for setting tooltip with plotly later to indicate hover text
(scores_noqc_plot <- scores_noqc %>%
  ggplot(aes(x = PC1, y = PC2, fill = treatment, text = glue("Sample: {sample_name},
                                                             Treatment: {treatment}"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 3, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5",
                                "QC" = "#F6A900")) +
  labs(x = glue("PC1: {PC1_percent_noqc}%"), 
       y = glue("PC2: {PC2_percent_noqc}%"), 
       title = "", tag = "A"))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-24-1.png)

Then make your scores plot ineractive so you can see who is who.

``` r
ggplotly(scores_noqc_plot, tooltip = "text")
```

<div class="plotly html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-5979381a1be227a30e79" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-5979381a1be227a30e79">{"x":{"data":[{"x":[-37.5794774514711,39.4713679058553],"y":[0,0],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,0],"y":[-15.2015826192488,22.0128709922432],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-34.0771662988653,-33.138053408564,-33.6248696491429,-33.7523148228631,-33.5013461934883,-33.6021156661427],"y":[-13.5100165459991,-12.9989625721584,-12.9891680136042,-13.150557149876,-13.1060804815934,-13.0523885008632],"text":["Sample: FJ-1,<br />Treatment: FJ","Sample: FJ-3,<br />Treatment: FJ","Sample: FJ-4,<br />Treatment: FJ","Sample: FJ-6,<br />Treatment: FJ","Sample: FJ-5,<br />Treatment: FJ","Sample: FJ-2,<br />Treatment: FJ"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(242,142,173,1)","opacity":1,"size":11.3385826771654,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"FJ","legendgroup":"FJ","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[35.738837734527,28.2295723781038,24.494397922074,22.1974321674877,24.9265230921361],"y":[-5.63336941620024,-3.04835087660838,-1.76967330049588,-0.614980515734691,-3.50900415924024],"text":["Sample: MEF-1,<br />Treatment: MEF","Sample: MEF-2,<br />Treatment: MEF","Sample: MEF-4,<br />Treatment: MEF","Sample: MEF-5,<br />Treatment: MEF","Sample: MEF-3,<br />Treatment: MEF"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(168,216,185,1)","opacity":1,"size":11.3385826771654,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF","legendgroup":"MEF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[35.9690567532495,34.2830511369859,32.8931534274563,33.7696560626205,34.6701171723752],"y":[-2.07747780130036,-2.58384517143169,0.876839245411983,0.0375772398898027,-0.630891789951512],"text":["Sample: MEF+SS-6,<br />Treatment: MEF+SS","Sample: MEF+SS-4,<br />Treatment: MEF+SS","Sample: MEF+SS-1,<br />Treatment: MEF+SS","Sample: MEF+SS-2,<br />Treatment: MEF+SS","Sample: MEF+SS-3,<br />Treatment: MEF+SS"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(127,179,213,1)","opacity":1,"size":11.3385826771654,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF+SS","legendgroup":"MEF+SS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-21.009537794891,-21.3034134180221,-21.1556269134966,-20.4291030331137,-21.5782506484264],"y":[18.9791697061589,19.562067890038,20.3213049189936,19.3480285658916,19.5497787286737],"text":["Sample: MT-2,<br />Treatment: MT","Sample: MT-5,<br />Treatment: MT","Sample: MT-4,<br />Treatment: MT","Sample: MT-6,<br />Treatment: MT","Sample: MT-3,<br />Treatment: MT"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(199,124,219,1)","opacity":1,"size":11.3385826771654,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MT","legendgroup":"MT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":40.1826484018265,"l":43.1050228310502},"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-37.5794774514711,39.4713679058553],"tickmode":"array","ticktext":["-20","0","20"],"tickvals":[-20,0,20],"categoryorder":"array","categoryarray":["-20","0","20"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"PC1: 73.4%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-15.2015826192488,22.0128709922432],"tickmode":"array","ticktext":["-10","0","10","20"],"tickvals":[-10,0,10,20],"categoryorder":"array","categoryarray":["-10","0","10","20"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"PC2: 11.9%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"title":{"text":"treatment","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"ceb43ad2f19":{"yintercept":{},"type":"scatter"},"ceb427a86218":{"xintercept":{}},"ceb445434c0":{"x":{},"y":{},"fill":{},"text":{}}},"cur_data":"ceb43ad2f19","visdat":{"ceb43ad2f19":["function (y) ","x"],"ceb427a86218":["function (y) ","x"],"ceb445434c0":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

``` r
ggsave(plot = scores_noqc_plot,
       filename = "pca_noqc_final_neg.jpeg", width = 140, height = 120, units = "mm")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-27-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

Top 10 features that influence PCA separation.

``` r
top10_features <- pca_noqc$rotation[, 1] %>%
  abs() %>%
  sort(decreasing = TRUE) %>%
  head(10) %>%
  names()

top10_scores <- pca_noqc$rotation[top10_features, 1]

# Creating a dataframe with top10 features and their scores
top10 <- data.frame(
  Feature = top10_features,
  Score = top10_scores)

# Display the dataframe
top10
```

                            Feature       Score
    431.0651_2.1848 431.0651_2.1848  0.11904532
    329.1605_3.9789 329.1605_3.9789  0.11300498
    159.1027_3.4799 159.1027_3.4799  0.11079952
    131.0714_2.3515 131.0714_2.3515  0.10703123
    517.1707_3.291   517.1707_3.291  0.10591679
    435.2232_3.4991 435.2232_3.4991  0.10429090
    133.0414_0.323   133.0414_0.323 -0.09680197
    323.033_2.7877   323.033_2.7877  0.09652671
    347.1343_2.5577 347.1343_2.5577  0.09168058
    521.2014_2.5235 521.2014_2.5235  0.08269297

## Add annotations

``` r
feature_table_long_log2 <- left_join(feature_table_long_log2, annotation, by = "mz_rt")
  
unique(feature_table_long_log2$Compound)
```

     [1] NA                                      
     [2] "2,4,6-Tryhydroxybenzaldehyde"          
     [3] "3,4-Dihydroxybenzoic acid"             
     [4] "p-Coumaric acid"                       
     [5] "2,4,6-Trihydroxybenzoic acid"          
     [6] "Gallic acid"                           
     [7] "Caffeic acid"                          
     [8] "Sorbitol"                              
     [9] "Ferulic acid"                          
    [10] "Pantothenic acid"                      
    [11] "Sinapic acid"                          
    [12] "Glucose"                               
    [13] "Phloretin"                             
    [14] "Epicatechin"                           
    [15] "Ellagic acid"                          
    [16] "Coumaroyl quinic acid"                 
    [17] "Sucrose"                               
    [18] "3-Caffeoylquinic acid"                 
    [19] "5-Caffeoylquinic acid"                 
    [20] "Feruloylquinic acid"                   
    [21] "Phlorhizin"                            
    [22] "Quercitrin"                            
    [23] "Quercetin-3-glucoside"                 
    [24] "Procyanidin B2"                        
    [25] "Naringin"                              
    [26] "Kampferol -3-O-glucosyl(1-2)rhamnoside"
    [27] "Eriocitrin"                            
    [28] "Rutin"                                 
    [29] "Procyanidin C1"                        

``` r
feature_table_long_log2 <- feature_table_long_log2 %>%
  mutate(id = if_else(is.na(Compound), mz_rt, Compound))
```

``` r
write_csv(feature_table_long_log2,
          file = "feature_table_long_log2__annotated_neg.csv")
```

## Univariate Testing

I am going to include some sample univariate testing for comparisons
between two groups.

### Fresh x MEF

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across fresh juice and MEF
freshjuice_MEF_novariance <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
freshjuice_MEF_toremove <- freshjuice_MEF_novariance$id

# investigating one of the features that has zero variance  
feature_table_long_log2 %>%
  filter(mz_rt == "116.0352_0.9495") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(treatment) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund))  
```

    # A tibble: 5 × 3
      treatment  stdev  mean
      <fct>      <dbl> <dbl>
    1 FJ        0       14.7
    2 MEF       0.0567  18.1
    3 MEF+SS    0.127   17.9
    4 MT        0       14.7
    5 QC        0.201   16.0

``` r
# run series of t-tests
freshjuice_MEF <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF") %>%
  filter(!mz_rt %in% freshjuice_MEF_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
freshjuice_MEF_bonferroni <- freshjuice_MEF %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_MEF_bonferroni)
```

    [1] 1445

#### Volcano plot

Wrangling

``` r
freshjuice_MEF_FC <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(freshjuice_minus_MEF_log2 = `FJ`-MEF)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
freshjuice_MEF_FC_pval <- right_join(freshjuice_MEF_FC, freshjuice_MEF, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_freshjuice <- freshjuice_MEF_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_MEF_log2 >= 1)

higher_MEF <- freshjuice_MEF_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_MEF_log2 <= -1)

(freshjuice_MEF_volcano <- freshjuice_MEF_FC_pval %>%
  ggplot(aes(x = freshjuice_minus_MEF_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_freshjuice, 
             aes(x = freshjuice_minus_MEF_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_MEF, 
             aes(x = freshjuice_minus_MEF_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "FJ vs MEF", tag = "B", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-33-1.png)

Make interactive.

``` r
ggplotly(freshjuice_MEF_volcano, tooltip = "text")
```

``` r
freshjuice_MEF_volcano <- freshjuice_MEF_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_freshjuice)} features higher in FJ"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_MEF)} features higher in MEF"), size = 2)
```

Saving higher_freshjuice and higher_MEF files

``` r
write_csv(higher_freshjuice,
          "neg_higher_freshjuice_versus_MEF.csv")

write_csv(higher_MEF,
          "neg_higher_MEF_versus_freshjuice.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
freshjuice_MEF_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF") %>%
  filter(!mz_rt %in% freshjuice_MEF_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
freshjuice_MEF_nonparametric_sig <- freshjuice_MEF_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_MEF_nonparametric_sig)
```

    [1] 1418

### Freshjuice x MEF_SS

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across fresh juice and MEF_SS
freshjuice_MEF_SS_novariance <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF+SS") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
freshjuice_MEF_SS_toremove <- freshjuice_MEF_SS_novariance$id

# investigating one of the features that has zero variance  
feature_table_long_log2 %>%
  filter(mz_rt == "116.0352_0.9495") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(treatment) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund))  
```

    # A tibble: 5 × 3
      treatment  stdev  mean
      <fct>      <dbl> <dbl>
    1 FJ        0       14.7
    2 MEF       0.0567  18.1
    3 MEF+SS    0.127   17.9
    4 MT        0       14.7
    5 QC        0.201   16.0

``` r
# run series of t-tests
freshjuice_MEF_SS <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF+SS") %>%
  filter(!mz_rt %in% freshjuice_MEF_SS_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "bonferroni") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
freshjuice_MEF_SS_bonferroni <- freshjuice_MEF_SS %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_MEF_SS_bonferroni)
```

    [1] 1566

#### Volcano plot

Wrangling

``` r
freshjuice_MEF_SS_FC <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF+SS") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(freshjuice_minus_MEF_SS_log2 = `FJ`-`MEF+SS`)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
freshjuice_MEF_SS_FC_pval <- right_join(freshjuice_MEF_SS_FC, freshjuice_MEF_SS, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_freshjuice <- freshjuice_MEF_SS_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_MEF_SS_log2 >= 1)

higher_MEF_SS <- freshjuice_MEF_SS_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_MEF_SS_log2 <= -1)

(freshjuice_MEF_SS_volcano <- freshjuice_MEF_SS_FC_pval %>%
  ggplot(aes(x = freshjuice_minus_MEF_SS_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_freshjuice, 
             aes(x = freshjuice_minus_MEF_SS_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_MEF_SS, 
             aes(x = freshjuice_minus_MEF_SS_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "FJ vs MEF+SS", tag = "C", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-40-1.png)

Make interactive.

``` r
ggplotly(freshjuice_MEF_SS_volcano, tooltip = "text")
```

``` r
freshjuice_MEF_SS_volcano <- freshjuice_MEF_SS_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_freshjuice)} features higher in FJ"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_MEF_SS)} features higher in MEF+SS"), size = 2)
```

Saving higher_freshjuice and higher_MEF_SS files

``` r
write_csv(higher_freshjuice,
          "neg_higher_freshjuice_versus_MEF_SS.csv")

write_csv(higher_MEF_SS,
          "neg_higher_MEF_SS_versus_freshjuice.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
freshjuice_MEF_SS_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MEF+SS") %>%
  filter(!mz_rt %in% freshjuice_MEF_SS_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
freshjuice_MEF_SS_nonparametric_sig <- freshjuice_MEF_SS_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_MEF_SS_nonparametric_sig)
```

    [1] 1534

### MEF x MEF_SS

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across MEF and MEF_SS
MEF_MEF_SS_novariance <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MEF+SS") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
MEF_MEF_SS_toremove <- MEF_MEF_SS_novariance$id

# investigating one of the features that has zero variance  
feature_table_long_log2 %>%
  filter(mz_rt == "116.0352_0.9495") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(treatment) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund))  
```

    # A tibble: 5 × 3
      treatment  stdev  mean
      <fct>      <dbl> <dbl>
    1 FJ        0       14.7
    2 MEF       0.0567  18.1
    3 MEF+SS    0.127   17.9
    4 MT        0       14.7
    5 QC        0.201   16.0

``` r
# run series of t-tests
MEF_MEF_SS <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MEF+SS") %>%
  filter(!mz_rt %in% MEF_MEF_SS_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
MEF_MEF_SS_bonferroni <- MEF_MEF_SS %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_MEF_SS_bonferroni)
```

    [1] 726

#### Volcano plot

Wrangling

``` r
MEF_MEF_SS_FC <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MEF+SS") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(MEF_minus_MEF_SS_log2 = MEF-`MEF+SS`)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
MEF_MEF_SS_FC_pval <- right_join(MEF_MEF_SS_FC, MEF_MEF_SS, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_MEF <- MEF_MEF_SS_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_minus_MEF_SS_log2 >= 1)

higher_MEF_SS <- MEF_MEF_SS_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_minus_MEF_SS_log2 <= -1)

(MEF_MEF_SS_volcano <- MEF_MEF_SS_FC_pval %>%
  ggplot(aes(x = MEF_minus_MEF_SS_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_MEF, 
             aes(x = MEF_minus_MEF_SS_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_MEF_SS, 
             aes(x = MEF_minus_MEF_SS_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +  
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "MEF vs MEF+SS", tag = "D", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-47-1.png)

Make interactive.

``` r
ggplotly(MEF_MEF_SS_volcano, tooltip = "text")
```

``` r
MEF_MEF_SS_volcano <- MEF_MEF_SS_volcano <- MEF_MEF_SS_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_MEF)} features higher in MEF"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_MEF_SS)} features higher in MEF+SS"), size = 2)
```

Saving higher_MEF and higher_MEF_SS files

``` r
write_csv(higher_MEF,
          "neg_higher_MEF_versus_MEF_SS.csv")

write_csv(higher_MEF_SS,
          "neg_higher_MEF_SS_versus_MEF.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
MEF_MEF_SS_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MEF+SS") %>%
  filter(!mz_rt %in% MEF_MEF_SS_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
MEF_MEF_SS_nonparametric_sig <- MEF_MEF_SS_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_MEF_SS_nonparametric_sig)
```

    [1] 593

### MEF x 40TT_10M

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across MEF and 40TT_10M
MEF_40TT_10M_novariance <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MT") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
MEF_40TT_10M_toremove <- MEF_40TT_10M_novariance$id


# run series of t-tests
MEF_40TT_10M <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MT") %>%
  filter(!mz_rt %in% MEF_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
MEF_40TT_10M_bonferroni <- MEF_40TT_10M %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_40TT_10M_bonferroni)
```

    [1] 1137

#### Volcano plot

Wrangling

``` r
MEF_40TT_10M_FC <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MT") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(MEF_minus_40TT_10M_log2 = `MEF`-`MT`)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
MEF_40TT_10M_FC_pval <- right_join(MEF_40TT_10M_FC, MEF_40TT_10M, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_MEF <- MEF_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_minus_40TT_10M_log2 >= 1)

higher_40TT_10M <- MEF_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_minus_40TT_10M_log2 <= -1)

(MEF_40TT_10M_volcano <- MEF_40TT_10M_FC_pval %>%
  ggplot(aes(x = MEF_minus_40TT_10M_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_MEF, 
             aes(x = MEF_minus_40TT_10M_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_40TT_10M, 
             aes(x = MEF_minus_40TT_10M_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +  
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "MEF vs MT", tag = "E", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-54-1.png)

Make interactive.

``` r
ggplotly(MEF_40TT_10M_volcano, tooltip = "text")
```

``` r
MEF_40TT_10M_volcano <- MEF_40TT_10M_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_MEF)} features higher in MEF"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_40TT_10M)} features higher in MT"), size = 2)
```

Saving higher_MEF and higher_40TT_10M files

``` r
write_csv(higher_MEF,
          "neg_higher_MEF_versus_40TT_10M.csv")

write_csv(higher_40TT_10M,
          "neg_higher_40TT_10M_versus_MEF.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
MEF_40TT_10M_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "MEF" | treatment == "MT") %>%
  filter(!mz_rt %in% MEF_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
MEF_40TT_10M_nonparametric_sig <- MEF_40TT_10M_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_40TT_10M_nonparametric_sig)
```

    [1] 1108

### MEF_SS x 40TT_10M

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across MEF_SS and 40TT_10M
MEF_SS_40TT_10M_novariance <- feature_table_long_log2 %>%
  filter(treatment == "MEF+SS" | treatment == "MT") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
MEF_SS_40TT_10M_toremove <- MEF_SS_40TT_10M_novariance$id


# run series of t-tests
MEF_SS_40TT_10M <- feature_table_long_log2 %>%
  filter(treatment == "MEF+SS" | treatment == "MT") %>%
  filter(!mz_rt %in% MEF_SS_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
MEF_SS_40TT_10M_bonferroni <- MEF_SS_40TT_10M %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_SS_40TT_10M_bonferroni)
```

    [1] 1341

#### Volcano plot

Wrangling

``` r
MEF_SS_40TT_10M_FC <- feature_table_long_log2 %>%
  filter(treatment == "MEF+SS" | treatment == "MT") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(MEF_SS_minus_40TT_10M_log2 = `MEF+SS`-`MT`)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
MEF_SS_40TT_10M_FC_pval <- right_join(MEF_SS_40TT_10M_FC, MEF_SS_40TT_10M, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_MEF_SS <- MEF_SS_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_SS_minus_40TT_10M_log2 >= 1)

higher_40TT_10M <- MEF_SS_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & MEF_SS_minus_40TT_10M_log2 <= -1)

(MEF_SS_40TT_10M_volcano <- MEF_SS_40TT_10M_FC_pval %>%
  ggplot(aes(x = MEF_SS_minus_40TT_10M_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_MEF_SS, 
             aes(x = MEF_SS_minus_40TT_10M_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_40TT_10M, 
             aes(x = MEF_SS_minus_40TT_10M_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +  
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "MEF+SS vs MT", tag = "F", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-61-1.png)

Make interactive.

``` r
ggplotly(MEF_SS_40TT_10M_volcano, tooltip = "text")
```

``` r
MEF_SS_40TT_10M_volcano <- MEF_SS_40TT_10M_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_MEF_SS)} features higher in MEF+SS"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_40TT_10M)} features higher in MT"), size = 2)
```

Saving higher_MEF_SS and higher_40TT_10M files

``` r
write_csv(higher_MEF_SS,
          "neg_higher_MEF_SS_versus_40TT_10M.csv")

write_csv(higher_40TT_10M,
          "neg_higher_40TT_10M_versus_MEF_SS.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
MEF_SS_40TT_10M_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "MEF+SS" | treatment == "MT") %>%
  filter(!mz_rt %in% MEF_SS_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
MEF_SS_40TT_10M_nonparametric_sig <- MEF_SS_40TT_10M_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(MEF_SS_40TT_10M_nonparametric_sig)
```

    [1] 1292

### Freshjuice x 40TT_10M

#### Parametric

Getting an error that some data is essential constant, checked for which
features have 0 SD across the treatments for comparison.

``` r
# find which features have no variance across freshjuice and 40TT_10M
freshjuice_40TT_10M_novariance <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MT") %>%
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# create a vector of those features with no variance
freshjuice_40TT_10M_toremove <- freshjuice_40TT_10M_novariance$id


# run series of t-tests
freshjuice_40TT_10M <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MT") %>%
  filter(!mz_rt %in% freshjuice_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  t_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # bonferroni or Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance()

# extract out only the significantly different features
freshjuice_40TT_10M_bonferroni <- freshjuice_40TT_10M %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_40TT_10M_bonferroni)
```

    [1] 827

#### Volcano plot

Wrangling

``` r
freshjuice_40TT_10M_FC <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MT") %>%
  group_by(treatment, id) %>%
  summarize(mean = mean(rel_abund)) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(freshjuice_minus_40TT_10M_log2 = `FJ`-`MT`)
```

    `summarise()` has grouped output by 'treatment'. You can override using the
    `.groups` argument.

``` r
freshjuice_40TT_10M_FC_pval <- right_join(freshjuice_40TT_10M_FC, freshjuice_40TT_10M, by = "id") %>%
  mutate(neglog10p = -log10(p))
```

Looking at features that are at least 2 fold change between groups and
significantly different at p\<0.05.

``` r
higher_freshjuice <- freshjuice_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_40TT_10M_log2 >= 1)

higher_40TT_10M <- freshjuice_40TT_10M_FC_pval %>%
  filter(neglog10p >= 1.3 & freshjuice_minus_40TT_10M_log2 <= -1)

(freshjuice_40TT_10M_volcano <- freshjuice_40TT_10M_FC_pval %>%
  ggplot(aes(x = freshjuice_minus_40TT_10M_log2, y = neglog10p, text = id)) +
  geom_point() +
  geom_point(data = higher_freshjuice, 
             aes(x = freshjuice_minus_40TT_10M_log2, y = neglog10p),
             color = "red") +
  geom_point(data = higher_40TT_10M, 
             aes(x = freshjuice_minus_40TT_10M_log2, y = neglog10p),
             color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(-8,8)) +
  theme_minimal(base_size = 8) +
  labs(x = "-log2(fold change)", base_size = 8,
       y = "-log10(p-value)", base_size = 8,
       title = "FJ vs MT", tag = "A", base_size = 8))
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-I_files/figure-commonmark/unnamed-chunk-68-1.png)

Make interactive.

``` r
ggplotly(freshjuice_40TT_10M_volcano, tooltip = "text")
```

``` r
freshjuice_40TT_10M_volcano <- freshjuice_40TT_10M_volcano +
  annotate(geom = "text", x = 4, y = 14, label = glue("{nrow(higher_freshjuice)} features higher in FJ"), size = 2) +
  annotate(geom = "text", x = -4, y = 14, label = glue("{nrow(higher_40TT_10M)} features higher in MT"), size = 2)
```

Saving higher_freshjuice and higher_40TT_10M files

``` r
write_csv(higher_freshjuice,
          "neg_higher_freshjuice_versus_40TT_10M.csv")

write_csv(higher_40TT_10M,
          "neg_higher_40TT_10M_versus_freshjuice.csv")
```

#### Non-parametric

Data might not be normally distributed so I did a nonparametric test.

``` r
# run series of t-tests
freshjuice_40TT_10M_nonparametric <- feature_table_long_log2 %>%
  filter(treatment == "FJ" | treatment == "MT") %>%
  filter(!mz_rt %in% freshjuice_40TT_10M_toremove) %>% # remove features with 0 variance
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
  wilcox_test(rel_abund ~ treatment, 
         paired = FALSE, 
         detailed = TRUE, # gives you more detail in output
         p.adjust.method = "BH") %>% # Benjamini-Hochberg false discovery rate multiple testing correction
  add_significance() %>%
  arrange(p)

# extract out only the significantly different features
freshjuice_40TT_10M_nonparametric_sig <- freshjuice_40TT_10M_nonparametric %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(freshjuice_40TT_10M_nonparametric_sig)
```

    [1] 768

``` r
volcano_plot <- (freshjuice_40TT_10M_volcano + freshjuice_MEF_volcano) / (freshjuice_MEF_SS_volcano + MEF_MEF_SS_volcano) / ( MEF_40TT_10M_volcano + MEF_SS_40TT_10M_volcano)

ggsave(plot = volcano_plot,
       filename = "volcano_plot_final_neg.svg", width = 190, height = 190, units = "mm")
```

## Heatmap of all features

``` r
head(feature_table_wide_log2)
```

    # A tibble: 6 × 2,289
      sample_name treatment run_order sample_or_qc `100.0402_0.8975`
      <chr>       <fct>         <dbl> <fct>                    <dbl>
    1 MT-2        MT               16 Sample                    16.8
    2 FJ-1        FJ               24 Sample                    16.7
    3 FJ-3        FJ               64 Sample                    16.9
    4 MT-5        MT               40 Sample                    16.7
    5 MT-4        MT               26 Sample                    16.7
    6 FJ-4        FJ               38 Sample                    16.9
    # ℹ 2,284 more variables: `100.9439_0.3382` <dbl>, `101.0222_0.3207` <dbl>,
    #   `101.0238_0.4781` <dbl>, `101.031_0.4766` <dbl>, `103.003_0.4435` <dbl>,
    #   `103.0384_0.3623` <dbl>, `103.0396_0.3544` <dbl>, `103.0397_0.6197` <dbl>,
    #   `103.0399_0.8555` <dbl>, `103.0399_1.01` <dbl>, `107.0348_0.9384` <dbl>,
    #   `111.0082_0.5501` <dbl>, `111.0087_0.3485` <dbl>, `111.0087_2.2117` <dbl>,
    #   `111.0088_1.9584` <dbl>, `111.0088_1.6771` <dbl>, `111.0139_0.3509` <dbl>,
    #   `111.0148_0.4951` <dbl>, `111.0156_0.6661` <dbl>, …

``` r
unique(feature_table_wide_log2$treatment)
```

    [1] MT     FJ     MEF+SS MEF    QC    
    Levels: FJ MEF MEF+SS MT QC

``` r
# find which features have no variance across fresh juice and MEF
# having samples with no variance won't work for a heatmap
all_MEF_comparisons_novariance <- feature_table_long_log2 %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# what features have no variance?
all_MEF_comparisons_novariance$id
```

    [1] "419.0805_4.825"  "636.5316_6.7983" "708.4612_6.7979"

``` r
# wrangle
feature_table_MEF_for_heatmap <- feature_table_wide_log2 %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  select(-`419.0805_4.825`, -`636.5316_6.7983`, -`708.4612_6.7979`) %>% # remove features with no variance
  select(-run_order, -sample_or_qc) %>% # remove columns we don't need
  column_to_rownames(var = "sample_name") # move the sample names to rownames for plotting


heatmap <- pheatmap(feature_table_MEF_for_heatmap[,-1], # remove treatment column
         cluster_rows = TRUE, # HCA on samples
         cluster_cols = TRUE, # HCA on features 
         scale = "column",
         legend_breaks = c(-4, 0, 4),
         legend_labels = c('Low relative abundance', 'Medium relative abundance', 'High relative abundance'),
         border_color = F,
         show_colnames = FALSE, # scale for each feature
         filename = "heatmap.png")
```

## ANOVA across groups

``` r
anova <- feature_table_long_log2 %>%
  filter(!sample_or_qc == "QC") %>% # remove QCs
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
 # Remove grupos onde todos os valores de rel_abund são iguais
  filter(n_distinct(rel_abund) > 1) %>%
  anova_test(rel_abund ~ treatment)

# adjust pvalues for multiple testing
anova_padjusted <- p.adjust(anova$p, method = "BH") %>%
  as.data.frame() %>%
  rename(p_adj = 1)

anova_padj <- bind_cols(as.data.frame(anova), anova_padjusted) %>%
  rename(padj = 9)

# extract out only the significantly different features
anova_sig <- anova_padj %>%
  as.data.frame() %>%
  filter(p <= 0.05)

# how many features are significantly different between the groups?
nrow(anova_sig)
```

    [1] 1713

## Heatmap of features that are different between groups

``` r
feature_table_long_log2_anova <- inner_join(feature_table_long_log2, anova_sig, by = "mz_rt") %>%
   select(sample_name, treatment, sample_or_qc, mz_rt, rel_abund, id)
  
unique(feature_table_long_log2$id)
```

       [1] "100.0402_0.8975"                       
       [2] "100.9439_0.3382"                       
       [3] "101.0222_0.3207"                       
       [4] "101.0238_0.4781"                       
       [5] "101.031_0.4766"                        
       [6] "103.003_0.4435"                        
       [7] "103.0384_0.3623"                       
       [8] "103.0396_0.3544"                       
       [9] "103.0397_0.6197"                       
      [10] "103.0399_0.8555"                       
      [11] "103.0399_1.01"                         
      [12] "107.0348_0.9384"                       
      [13] "111.0082_0.5501"                       
      [14] "111.0087_0.3485"                       
      [15] "111.0087_2.2117"                       
      [16] "111.0088_1.9584"                       
      [17] "111.0088_1.6771"                       
      [18] "111.0139_0.3509"                       
      [19] "111.0148_0.4951"                       
      [20] "111.0156_0.6661"                       
      [21] "111.0157_1.6749"                       
      [22] "111.0163_1.5561"                       
      [23] "111.0332_0.6322"                       
      [24] "112.0116_0.3425"                       
      [25] "113.0131_1.6741"                       
      [26] "113.0221_0.3115"                       
      [27] "113.0313_0.4599"                       
      [28] "115.0031_0.3715"                       
      [29] "115.0035_0.7176"                       
      [30] "115.0036_0.902"                        
      [31] "115.0037_1.5429"                       
      [32] "115.0096_0.3808"                       
      [33] "115.0104_0.3883"                       
      [34] "115.0399_0.7668"                       
      [35] "115.04_1.0671"                         
      [36] "115.04_2.1988"                         
      [37] "115.0401_3.6659"                       
      [38] "116.0063_0.3786"                       
      [39] "116.0345_0.4908"                       
      [40] "116.0352_0.9495"                       
      [41] "117.0063_0.4101"                       
      [42] "117.0076_0.3818"                       
      [43] "117.0178_0.3615"                       
      [44] "117.019_0.3477"                        
      [45] "117.0193_0.78"                         
      [46] "117.0265_0.3557"                       
      [47] "117.0557_1.817"                        
      [48] "117.0557_1.9458"                       
      [49] "125.0224_0.3609"                       
      [50] "125.0243_0.7922"                       
      [51] "125.0244_1.5476"                       
      [52] "126.9352_0.592"                        
      [53] "129.0187_0.4854"                       
      [54] "129.0191_0.6312"                       
      [55] "129.0193_0.7456"                       
      [56] "129.0193_1.0397"                       
      [57] "129.0193_1.5556"                       
      [58] "129.0254_0.3095"                       
      [59] "129.0265_0.4811"                       
      [60] "129.0267_0.3671"                       
      [61] "129.0274_1.0365"                       
      [62] "129.0558_2.1945"                       
      [63] "130.0206_0.3173"                       
      [64] "130.0221_0.4793"                       
      [65] "130.0225_0.6229"                       
      [66] "130.0227_1.0291"                       
      [67] "130.0227_0.3054"                       
      [68] "130.9981_0.6289"                       
      [69] "131.0344_0.4711"                       
      [70] "131.035_1.7432"                        
      [71] "131.0434_1.7452"                       
      [72] "131.0714_2.3515"                       
      [73] "131.0714_2.5701"                       
      [74] "131.0714_3.0302"                       
      [75] "133.0133_0.3672"                       
      [76] "133.0142_0.7621"                       
      [77] "133.0143_0.9193"                       
      [78] "133.0143_1.5273"                       
      [79] "133.0147_0.3722"                       
      [80] "133.0209_0.4699"                       
      [81] "133.0223_0.7556"                       
      [82] "133.0224_0.9214"                       
      [83] "133.0406_0.4428"                       
      [84] "133.0414_0.323"                        
      [85] "133.0505_0.6894"                       
      [86] "133.0506_1.0344"                       
      [87] "133.1535_0.3813"                       
      [88] "134.0173_0.3814"                       
      [89] "134.0176_0.7554"                       
      [90] "134.0176_0.9175"                       
      [91] "134.025_0.3857"                        
      [92] "135.0181_0.3827"                       
      [93] "137.0245_1.6923"                       
      [94] "139.0401_3.8511"                       
      [95] "142.0332_1.7694"                       
      [96] "142.0509_0.8979"                       
      [97] "142.051_1.5761"                        
      [98] "143.0328_0.325"                        
      [99] "143.0714_3.2474"                       
     [100] "143.0714_2.6454"                       
     [101] "144.0666_1.1017"                       
     [102] "145.0125_0.3651"                       
     [103] "145.0136_0.5034"                       
     [104] "145.0142_1.6178"                       
     [105] "145.0506_1.5905"                       
     [106] "145.0871_3.5012"                       
     [107] "147.0295_0.6175"                       
     [108] "147.0299_1.0983"                       
     [109] "147.0299_1.5736"                       
     [110] "147.0306_0.8544"                       
     [111] "147.0378_0.9012"                       
     [112] "147.0452_3.0243"                       
     [113] "147.0582_0.8554"                       
     [114] "147.0662_1.016"                        
     [115] "147.0663_1.5601"                       
     [116] "147.1767_0.8534"                       
     [117] "148.0333_0.8554"                       
     [118] "148.0418_0.8537"                       
     [119] "149.0245_1.792"                        
     [120] "149.0343_0.8536"                       
     [121] "149.0608_6.1724"                       
     [122] "151.0401_1.899"                        
     [123] "151.0401_2.3989"                       
     [124] "153.0193_2.4502"                       
     [125] "2,4,6-Tryhydroxybenzaldehyde"          
     [126] "3,4-Dihydroxybenzoic acid"             
     [127] "153.0322_1.1114"                       
     [128] "153.0322_1.3243"                       
     [129] "154.998_0.4815"                        
     [130] "154.9986_1.5388"                       
     [131] "156.0666_1.8198"                       
     [132] "157.0506_2.0254"                       
     [133] "157.0619_1.6608"                       
     [134] "159.0293_0.507"                        
     [135] "159.0298_0.7712"                       
     [136] "159.0299_1.7093"                       
     [137] "159.0387_0.7664"                       
     [138] "159.1027_3.131"                        
     [139] "159.1027_3.4799"                       
     [140] "160.0332_0.7611"                       
     [141] "161.0244_2.2121"                       
     [142] "161.0435_0.3562"                       
     [143] "161.0456_1.724"                        
     [144] "161.0545_1.7257"                       
     [145] "161.0819_2.6781"                       
     [146] "161.0819_2.0134"                       
     [147] "163.04_2.0089"                         
     [148] "163.0401_3.2283"                       
     [149] "163.0401_2.2402"                       
     [150] "163.0401_2.8135"                       
     [151] "p-Coumaric acid"                       
     [152] "164.0718_1.6899"                       
     [153] "165.04_0.3064"                         
     [154] "165.0761_0.3314"                       
     [155] "167.035_2.993"                         
     [156] "167.035_1.8278"                        
     [157] "167.035_2.7633"                        
     [158] "167.035_2.284"                         
     [159] "168.0302_1.6238"                       
     [160] "2,4,6-Trihydroxybenzoic acid"          
     [161] "Gallic acid"                           
     [162] "170.9612_1.6471"                       
     [163] "171.0663_2.993"                        
     [164] "171.0663_2.2588"                       
     [165] "173.0084_0.3557"                       
     [166] "173.0091_2.2149"                       
     [167] "173.0092_1.3425"                       
     [168] "173.0183_0.6803"                       
     [169] "173.0184_1.3431"                       
     [170] "173.0396_0.3482"                       
     [171] "173.0453_1.1037"                       
     [172] "173.0456_2.2521"                       
     [173] "173.0456_2.4643"                       
     [174] "173.0816_0.7892"                       
     [175] "173.0819_1.7242"                       
     [176] "173.0819_3.0287"                       
     [177] "173.0821_2.0891"                       
     [178] "173.1677_0.4796"                       
     [179] "174.0123_0.6322"                       
     [180] "174.0126_1.6739"                       
     [181] "174.0405_0.48"                         
     [182] "174.9859_0.3067"                       
     [183] "175.013_0.4806"                        
     [184] "175.0133_0.6139"                       
     [185] "175.0235_0.3294"                       
     [186] "175.0246_0.8599"                       
     [187] "175.0247_1.1048"                       
     [188] "175.0249_2.7543"                       
     [189] "175.0612_1.7219"                       
     [190] "175.0613_2.1993"                       
     [191] "175.0704_2.1996"                       
     [192] "177.0193_2.3613"                       
     [193] "177.0403_0.6944"                       
     [194] "Caffeic acid"                          
     [195] "179.0556_0.4552"                       
     [196] "179.056_0.7452"                        
     [197] "179.0628_0.3125"                       
     [198] "179.0648_0.458"                        
     [199] "180.0664_0.9078"                       
     [200] "181.0506_1.1034"                       
     [201] "181.0506_2.3328"                       
     [202] "181.0506_3.2741"                       
     [203] "Sorbitol"                              
     [204] "181.0806_0.4684"                       
     [205] "181.0811_0.7455"                       
     [206] "182.0459_1.5739"                       
     [207] "182.0745_0.4637"                       
     [208] "183.1027_4.2919"                       
     [209] "187.0248_1.5882"                       
     [210] "187.0612_1.7221"                       
     [211] "187.0612_1.8998"                       
     [212] "187.1339_4.2471"                       
     [213] "187.1341_3.603"                        
     [214] "189.0023_0.3273"                       
     [215] "189.0404_1.9197"                       
     [216] "189.0404_1.6168"                       
     [217] "189.0405_1.0482"                       
     [218] "189.0768_2.0488"                       
     [219] "189.0769_1.8561"                       
     [220] "189.0769_2.6695"                       
     [221] "191.0204_0.3473"                       
     [222] "191.0198_1.0624"                       
     [223] "191.0199_0.9072"                       
     [224] "191.0199_1.5621"                       
     [225] "191.0211_0.635"                        
     [226] "191.0281_0.815"                        
     [227] "191.0286_1.5634"                       
     [228] "191.0294_2.2139"                       
     [229] "191.0508_1.5665"                       
     [230] "191.0516_0.4748"                       
     [231] "191.052_0.7027"                        
     [232] "191.0528_1.0744"                       
     [233] "191.0548_1.6482"                       
     [234] "191.0561_2.4461"                       
     [235] "191.0563_2.2115"                       
     [236] "191.0913_0.3602"                       
     [237] "191.1876_0.3528"                       
     [238] "191.1879_0.6192"                       
     [239] "191.2049_0.3434"                       
     [240] "191.2047_0.6127"                       
     [241] "191.2171_0.3395"                       
     [242] "191.2171_0.624"                        
     [243] "191.231_0.6226"                        
     [244] "191.2313_0.3434"                       
     [245] "191.2488_0.6146"                       
     [246] "191.249_0.3543"                        
     [247] "191.2671_0.6171"                       
     [248] "191.2673_0.3452"                       
     [249] "191.2814_0.6076"                       
     [250] "191.3052_0.6169"                       
     [251] "191.3784_0.6105"                       
     [252] "191.9375_0.622"                        
     [253] "191.9502_0.6199"                       
     [254] "191.9804_0.6297"                       
     [255] "191.9944_0.6163"                       
     [256] "192.0042_0.6331"                       
     [257] "192.0227_0.3434"                       
     [258] "192.0232_2.2137"                       
     [259] "192.0232_1.5643"                       
     [260] "192.0307_0.3454"                       
     [261] "192.0313_0.4604"                       
     [262] "192.0314_0.5879"                       
     [263] "192.0329_1.5635"                       
     [264] "192.0531_0.6121"                       
     [265] "192.0537_0.3341"                       
     [266] "192.055_0.6331"                        
     [267] "192.0953_0.6316"                       
     [268] "193.0238_0.3471"                       
     [269] "193.024_0.6314"                        
     [270] "193.0242_1.5641"                       
     [271] "193.0334_0.632"                        
     [272] "Ferulic acid"                          
     [273] "193.0507_2.1269"                       
     [274] "193.0536_0.6323"                       
     [275] "193.0608_2.2125"                       
     [276] "193.087_3.4512"                        
     [277] "194.0268_0.3516"                       
     [278] "194.027_0.6304"                        
     [279] "194.8758_2.9697"                       
     [280] "195.0299_2.9182"                       
     [281] "195.0506_0.3009"                       
     [282] "195.0509_0.7424"                       
     [283] "195.1391_2.4854"                       
     [284] "197.0087_1.1043"                       
     [285] "197.0089_2.8079"                       
     [286] "200.9717_0.7779"                       
     [287] "201.1133_4.099"                        
     [288] "202.0721_1.5324"                       
     [289] "202.9754_2.2028"                       
     [290] "202.9755_1.6658"                       
     [291] "203.0827_2.0022"                       
     [292] "203.0925_2.6794"                       
     [293] "205.0354_0.8751"                       
     [294] "205.0357_1.5424"                       
     [295] "205.0446_1.5426"                       
     [296] "205.0454_1.0986"                       
     [297] "205.0506_2.3943"                       
     [298] "206.0388_1.0995"                       
     [299] "206.0489_1.5437"                       
     [300] "206.0823_2.7149"                       
     [301] "206.9684_0.3619"                       
     [302] "207.0128_0.3103"                       
     [303] "207.0874_1.8041"                       
     [304] "209.065_0.3111"                        
     [305] "211.0006_2.6906"                       
     [306] "211.0248_2.9418"                       
     [307] "212.9717_1.6717"                       
     [308] "214.0721_1.8857"                       
     [309] "214.8889_2.566"                        
     [310] "215.0331_0.3023"                       
     [311] "215.0327_0.7679"                       
     [312] "215.0328_0.9772"                       
     [313] "215.0423_0.4603"                       
     [314] "215.0429_0.9902"                       
     [315] "215.043_1.1918"                        
     [316] "216.0361_0.7606"                       
     [317] "216.0362_0.9939"                       
     [318] "216.0877_2.1095"                       
     [319] "216.0878_1.8112"                       
     [320] "216.9663_3.7041"                       
     [321] "216.9812_1.82"                         
     [322] "217.0298_0.3003"                       
     [323] "217.0298_0.978"                        
     [324] "217.0299_0.7667"                       
     [325] "217.0482_0.7645"                       
     [326] "217.0483_0.9711"                       
     [327] "217.0484_1.4877"                       
     [328] "217.0717_2.023"                        
     [329] "218.0653_0.3061"                       
     [330] "Pantothenic acid"                      
     [331] "218.1136_1.8203"                       
     [332] "219.0497_1.0737"                       
     [333] "219.0511_2.0026"                       
     [334] "221.0293_0.3083"                       
     [335] "221.0385_0.3081"                       
     [336] "221.0642_0.3071"                       
     [337] "223.0449_0.3069"                       
     [338] "Sinapic acid"                          
     [339] "223.0612_2.2279"                       
     [340] "Glucose"                               
     [341] "225.0609_0.4682"                       
     [342] "225.0881_1.8024"                       
     [343] "227.0561_0.9622"                       
     [344] "227.069_1.715"                         
     [345] "229.0452_0.3151"                       
     [346] "229.0982_2.3473"                       
     [347] "229.1194_2.0098"                       
     [348] "230.0298_0.474"                        
     [349] "230.067_1.6814"                        
     [350] "231.0298_3.5987"                       
     [351] "233.0665_1.1029"                       
     [352] "235.0611_2.1476"                       
     [353] "235.0612_2.9052"                       
     [354] "236.9258_0.8515"                       
     [355] "236.9258_1.1027"                       
     [356] "236.9363_2.7769"                       
     [357] "239.0765_0.4963"                       
     [358] "241.0846_2.0658"                       
     [359] "243.0623_0.9515"                       
     [360] "243.1237_3.3864"                       
     [361] "244.0659_0.3019"                       
     [362] "244.9281_0.8588"                       
     [363] "244.998_1.7264"                        
     [364] "245.0302_1.681"                        
     [365] "245.0666_1.7441"                       
     [366] "245.0932_2.8687"                       
     [367] "246.0982_2.2753"                       
     [368] "246.0982_2.4891"                       
     [369] "246.9258_0.8574"                       
     [370] "246.9918_1.9463"                       
     [371] "247.0555_0.3146"                       
     [372] "247.0823_1.8226"                       
     [373] "248.0404_0.4735"                       
     [374] "248.1292_3.9102"                       
     [375] "249.025_0.6922"                        
     [376] "249.0615_3.4511"                       
     [377] "249.0979_1.638"                        
     [378] "249.0979_1.8046"                       
     [379] "249.1244_1.8566"                       
     [380] "251.1135_1.9333"                       
     [381] "251.1136_1.5886"                       
     [382] "252.0877_2.0782"                       
     [383] "253.0718_2.7214"                       
     [384] "253.0928_0.9282"                       
     [385] "253.1038_0.929"                        
     [386] "254.0962_0.9287"                       
     [387] "256.9613_1.5337"                       
     [388] "256.9613_1.672"                        
     [389] "256.9613_0.9735"                       
     [390] "259.1299_2.0058"                       
     [391] "260.0451_0.9456"                       
     [392] "261.0615_1.6366"                       
     [393] "261.088_2.3477"                        
     [394] "261.1345_2.1181"                       
     [395] "262.0558_0.3121"                       
     [396] "262.0673_0.3153"                       
     [397] "262.9043_0.4754"                       
     [398] "262.9051_1.6744"                       
     [399] "262.9051_1.5311"                       
     [400] "263.1288_2.8921"                       
     [401] "263.1288_3.5393"                       
     [402] "265.0928_1.1048"                       
     [403] "265.0929_1.5914"                       
     [404] "265.1039_1.5916"                       
     [405] "266.0961_1.5919"                       
     [406] "266.9078_1.7934"                       
     [407] "267.1085_1.6681"                       
     [408] "268.0825_1.6751"                       
     [409] "270.9057_0.4776"                       
     [410] "270.9758_0.4779"                       
     [411] "271.1663_2.5384"                       
     [412] "272.9339_1.5322"                       
     [413] "273.0244_1.5343"                       
     [414] "273.025_1.112"                         
     [415] "273.0355_0.5067"                       
     [416] "Phloretin"                             
     [417] "273.0881_3.1766"                       
     [418] "274.0246_0.423"                        
     [419] "274.972_0.7363"                        
     [420] "274.9721_0.971"                        
     [421] "274.9721_1.5738"                       
     [422] "275.021_0.4243"                        
     [423] "275.032_0.4157"                        
     [424] "275.1136_2.7039"                       
     [425] "275.1288_5.5664"                       
     [426] "275.1288_4.603"                        
     [427] "276.0175_0.9465"                       
     [428] "276.0245_0.4194"                       
     [429] "276.9287_0.7446"                       
     [430] "277.0564_1.635"                        
     [431] "277.1445_5.4452"                       
     [432] "278.0879_0.704"                        
     [433] "279.1602_4.3139"                       
     [434] "280.9155_0.7341"                       
     [435] "281.0698_2.2026"                       
     [436] "281.1241_2.1077"                       
     [437] "281.1393_2.4994"                       
     [438] "281.1759_4.1679"                       
     [439] "282.0844_1.5791"                       
     [440] "283.0823_2.3919"                       
     [441] "283.155_3.1958"                        
     [442] "284.0008_1.6926"                       
     [443] "284.1584_3.219"                        
     [444] "285.0402_3.5742"                       
     [445] "285.1109_2.5213"                       
     [446] "285.1343_5.2638"                       
     [447] "285.1706_2.9951"                       
     [448] "287.0408_1.62"                         
     [449] "287.0771_1.109"                        
     [450] "287.0902_1.59"                         
     [451] "287.1018_1.59"                         
     [452] "288.9179_0.7606"                       
     [453] "288.965_0.845"                         
     [454] "289.017_0.3931"                        
     [455] "289.0318_0.3153"                       
     [456] "289.0352_2.4576"                       
     [457] "289.0564_1.0931"                       
     [458] "Epicatechin"                           
     [459] "289.0823_2.435"                        
     [460] "289.0829_2.8729"                       
     [461] "289.0874_1.5898"                       
     [462] "289.1292_3.1598"                       
     [463] "289.1444_5.0488"                       
     [464] "289.1445_6.1724"                       
     [465] "289.156_5.0489"                        
     [466] "289.156_6.1727"                        
     [467] "290.0864_2.4293"                       
     [468] "290.9156_0.7579"                       
     [469] "291.0035_1.8418"                       
     [470] "291.16_4.7973"                         
     [471] "291.1965_5.9545"                       
     [472] "293.1241_1.9339"                       
     [473] "293.1241_2.2004"                       
     [474] "293.2121_6.2337"                       
     [475] "293.2121_5.7726"                       
     [476] "294.9312_1.6664"                       
     [477] "295.0663_0.3141"                       
     [478] "295.0669_0.7184"                       
     [479] "295.0776_0.3153"                       
     [480] "295.0823_2.5206"                       
     [481] "295.0826_2.1481"                       
     [482] "295.1033_1.1032"                       
     [483] "295.1034_1.5644"                       
     [484] "295.1035_1.8045"                       
     [485] "295.115_1.8041"                        
     [486] "295.1187_2.5911"                       
     [487] "295.1398_2.521"                        
     [488] "295.2277_6.0627"                       
     [489] "296.0697_0.4847"                       
     [490] "297.1192_1.5898"                       
     [491] "297.1306_1.5888"                       
     [492] "299.0771_2.1362"                       
     [493] "299.0772_1.6933"                       
     [494] "299.0892_1.6934"                       
     [495] "299.0902_2.2938"                       
     [496] "299.124_1.5891"                        
     [497] "299.1499_2.4852"                       
     [498] "299.973_1.6927"                        
     [499] "Ellagic acid"                          
     [500] "301.0717_3.0615"                       
     [501] "301.0927_1.6816"                       
     [502] "301.1655_2.9075"                       
     [503] "303.006_0.8012"                        
     [504] "303.0509_2.452"                        
     [505] "303.1448_3.5129"                       
     [506] "303.2176_3.5466"                       
     [507] "304.9908_0.3931"                       
     [508] "305.0666_2.2424"                       
     [509] "305.0667_1.7959"                       
     [510] "305.1758_6.1936"                       
     [511] "307.012_1.8052"                        
     [512] "307.0305_0.7002"                       
     [513] "307.0305_1.1033"                       
     [514] "307.0823_2.5985"                       
     [515] "307.1033_2.1358"                       
     [516] "307.1398_2.6147"                       
     [517] "307.1398_2.1875"                       
     [518] "307.1681_3.1868"                       
     [519] "307.1915_5.0995"                       
     [520] "307.1916_4.8151"                       
     [521] "308.0985_0.6981"                       
     [522] "308.9316_0.47"                         
     [523] "308.9317_0.354"                        
     [524] "309.1191_2.2928"                       
     [525] "309.1554_3.0407"                       
     [526] "310.1506_0.9387"                       
     [527] "311.1136_2.2211"                       
     [528] "311.1346_1.7272"                       
     [529] "312.0949_1.5583"                       
     [530] "313.0564_1.9308"                       
     [531] "313.1132_0.4684"                       
     [532] "313.1213_4.644"                        
     [533] "313.1292_2.3057"                       
     [534] "313.1444_6.1987"                       
     [535] "313.1655_3.9332"                       
     [536] "314.088_1.7952"                        
     [537] "315.0143_3.1437"                       
     [538] "315.072_2.1296"                        
     [539] "315.0873_2.2346"                       
     [540] "315.1084_1.8923"                       
     [541] "315.1811_3.8553"                       
     [542] "317.0301_2.5345"                       
     [543] "317.0488_0.8569"                       
     [544] "317.049_1.1014"                        
     [545] "317.0846_1.8039"                       
     [546] "317.1756_6.8643"                       
     [547] "317.1757_5.7877"                       
     [548] "319.0249_0.8565"                       
     [549] "319.0459_2.0537"                       
     [550] "319.0589_2.3928"                       
     [551] "319.0823_2.4332"                       
     [552] "320.092_1.0428"                        
     [553] "320.9542_0.4081"                       
     [554] "321.0383_1.6262"                       
     [555] "321.0627_0.3586"                       
     [556] "321.0825_0.7601"                       
     [557] "321.0826_1.7945"                       
     [558] "323.0277_0.4883"                       
     [559] "323.033_2.7877"                        
     [560] "323.0538_0.8087"                       
     [561] "323.0771_2.4103"                       
     [562] "323.0982_1.0971"                       
     [563] "323.0983_1.7702"                       
     [564] "323.0983_1.5656"                       
     [565] "323.1105_1.7696"                       
     [566] "323.1346_2.5602"                       
     [567] "323.1347_2.0298"                       
     [568] "323.1347_1.7622"                       
     [569] "323.1711_1.7244"                       
     [570] "324.1017_1.5636"                       
     [571] "325.0485_2.4276"                       
     [572] "325.0927_2.0084"                       
     [573] "325.1137_1.6489"                       
     [574] "325.1138_1.1039"                       
     [575] "325.1502_2.5082"                       
     [576] "325.1503_1.9666"                       
     [577] "327.072_2.1462"                        
     [578] "327.072_2.4262"                        
     [579] "327.092_0.309"                         
     [580] "327.0983_2.2568"                       
     [581] "327.1055_1.1019"                       
     [582] "327.1083_2.4026"                       
     [583] "327.1085_2.1272"                       
     [584] "327.1207_2.1272"                       
     [585] "327.2175_3.947"                        
     [586] "329.0116_0.4606"                       
     [587] "329.0877_2.5189"                       
     [588] "329.0877_1.6009"                       
     [589] "329.0877_1.1309"                       
     [590] "329.0877_1.7965"                       
     [591] "329.0999_1.7986"                       
     [592] "329.124_2.0908"                        
     [593] "329.1604_3.6089"                       
     [594] "329.1605_3.9789"                       
     [595] "330.0906_2.5183"                       
     [596] "330.091_1.5996"                        
     [597] "330.091_1.7992"                        
     [598] "330.0911_2.2818"                       
     [599] "331.028_0.7749"                        
     [600] "331.0589_2.1478"                       
     [601] "331.0668_1.1048"                       
     [602] "331.067_1.5848"                        
     [603] "331.103_1.825"                         
     [604] "331.1913_6.4054"                       
     [605] "331.1914_6.1596"                       
     [606] "331.1914_6.0127"                       
     [607] "333.0216_0.8563"                       
     [608] "333.0263_0.4269"                       
     [609] "333.0353_0.8526"                       
     [610] "333.0614_3.0459"                       
     [611] "333.1186_2.2678"                       
     [612] "333.2069_6.0806"                       
     [613] "333.2069_6.8668"                       
     [614] "334.025_0.8501"                        
     [615] "334.0845_0.8243"                       
     [616] "335.0538_1.693"                        
     [617] "335.077_2.4918"                        
     [618] "335.0791_0.8448"                       
     [619] "335.1345_2.8142"                       
     [620] "335.1345_2.4943"                       
     [621] "335.1346_1.0988"                       
     [622] "335.1346_1.7474"                       
     [623] "335.2226_5.8376"                       
     [624] "337.051_1.6928"                        
     [625] "337.0769_0.3274"                       
     [626] "337.077_0.5879"                        
     [627] "337.0774_1.0941"                       
     [628] "Coumaroyl quinic acid"                 
     [629] "337.104_2.4639"                        
     [630] "337.1048_2.6265"                       
     [631] "337.1139_1.7046"                       
     [632] "337.1141_2.1289"                       
     [633] "337.1261_2.1292"                       
     [634] "337.1503_2.9328"                       
     [635] "338.0962_2.1747"                       
     [636] "338.1084_2.4631"                       
     [637] "339.0552_0.3058"                       
     [638] "339.1084_2.5226"                       
     [639] "339.1191_2.1292"                       
     [640] "339.1295_1.978"                        
     [641] "339.1295_1.6256"                       
     [642] "339.1659_2.2179"                       
     [643] "340.1037_2.5517"                       
     [644] "340.1329_1.9666"                       
     [645] "341.0877_2.5759"                       
     [646] "341.0898_1.7182"                       
     [647] "341.1012_2.1788"                       
     [648] "Sucrose"                               
     [649] "341.1087_1.8218"                       
     [650] "341.1088_0.947"                        
     [651] "341.1239_2.4411"                       
     [652] "341.1452_1.9403"                       
     [653] "343.0669_1.8715"                       
     [654] "343.1032_2.0912"                       
     [655] "343.1032_2.1825"                       
     [656] "343.1033_1.7607"                       
     [657] "343.1156_1.9846"                       
     [658] "343.1397_2.406"                        
     [659] "343.1527_2.9264"                       
     [660] "345.1189_1.8565"                       
     [661] "345.1341_0.9434"                       
     [662] "345.1343_3.0847"                       
     [663] "345.1504_2.9264"                       
     [664] "345.1552_3.2788"                       
     [665] "345.1553_2.4711"                       
     [666] "345.1553_2.7726"                       
     [667] "345.2069_6.4762"                       
     [668] "346.0557_1.5574"                       
     [669] "346.078_1.5973"                        
     [670] "346.1587_2.7732"                       
     [671] "347.0222_0.4236"                       
     [672] "347.0981_1.6462"                       
     [673] "347.1343_2.5577"                       
     [674] "347.171_2.897"                         
     [675] "347.9783_0.8556"                       
     [676] "347.9784_1.1022"                       
     [677] "348.0757_1.7194"                       
     [678] "348.986_0.8567"                        
     [679] "349.0473_0.7077"                       
     [680] "349.0695_1.8993"                       
     [681] "349.0927_3.6153"                       
     [682] "349.1059_2.3061"                       
     [683] "349.1139_2.3756"                       
     [684] "349.1866_2.8739"                       
     [685] "349.1867_3.0958"                       
     [686] "349.2383_6.3293"                       
     [687] "349.9891_0.8558"                       
     [688] "351.0719_2.2176"                       
     [689] "351.072_3.2314"                        
     [690] "351.1295_1.9307"                       
     [691] "351.1301_2.4794"                       
     [692] "351.1407_2.479"                        
     [693] "351.1658_3.3877"                       
     [694] "351.1659_2.621"                        
     [695] "351.1709_2.4796"                       
     [696] "352.1454_2.4792"                       
     [697] "353.0719_0.3201"                       
     [698] "353.0719_0.5132"                       
     [699] "353.0843_0.3213"                       
     [700] "353.0844_0.5118"                       
     [701] "3-Caffeoylquinic acid"                 
     [702] "5-Caffeoylquinic acid"                 
     [703] "353.0999_1.9741"                       
     [704] "353.1003_2.3986"                       
     [705] "353.1088_1.6169"                       
     [706] "353.1088_0.4625"                       
     [707] "353.1096_0.3171"                       
     [708] "353.1309_2.2124"                       
     [709] "353.1347_2.479"                        
     [710] "353.1452_1.858"                        
     [711] "353.1815_2.6272"                       
     [712] "353.1817_2.9262"                       
     [713] "353.1941_2.9261"                       
     [714] "354.0753_0.5111"                       
     [715] "354.1033_2.2239"                       
     [716] "354.1038_1.9726"                       
     [717] "354.9992_2.0043"                       
     [718] "355.0669_1.9649"                       
     [719] "355.0675_1.8484"                       
     [720] "355.0862_0.4994"                       
     [721] "355.0868_0.302"                        
     [722] "355.1006_0.3024"                       
     [723] "355.1032_2.1556"                       
     [724] "355.1033_2.9525"                       
     [725] "355.1034_2.3678"                       
     [726] "355.1051_1.8956"                       
     [727] "355.1158_2.3675"                       
     [728] "355.1162_2.1897"                       
     [729] "356.1032_2.2211"                       
     [730] "356.9808_0.8541"                       
     [731] "356.996_3.1349"                        
     [732] "357.1027_0.3063"                       
     [733] "357.1088_2.3681"                       
     [734] "357.1156_0.3054"                       
     [735] "357.1188_2.2314"                       
     [736] "357.1189_1.64"                         
     [737] "357.1376_4.06"                         
     [738] "358.0264_0.3463"                       
     [739] "358.0269_0.6905"                       
     [740] "358.106_0.3062"                        
     [741] "359.0299_0.6899"                       
     [742] "359.0982_1.9345"                       
     [743] "359.1343_2.7541"                       
     [744] "359.1345_2.0719"                       
     [745] "359.1346_2.3053"                       
     [746] "360.0233_0.6888"                       
     [747] "360.0427_1.0091"                       
     [748] "360.2754_4.8575"                       
     [749] "360.949_0.4608"                        
     [750] "361.0384_0.8245"                       
     [751] "361.0696_2.0073"                       
     [752] "361.1068_0.9416"                       
     [753] "361.1139_2.0013"                       
     [754] "361.1138_1.774"                        
     [755] "361.1502_2.3265"                       
     [756] "361.1502_2.8555"                       
     [757] "361.1503_2.1307"                       
     [758] "361.1503_1.9841"                       
     [759] "362.291_5.0366"                        
     [760] "362.9911_0.7452"                       
     [761] "362.9962_0.4067"                       
     [762] "363.0148_0.8219"                       
     [763] "363.0851_2.1274"                       
     [764] "363.0931_1.7254"                       
     [765] "363.1294_2.6905"                       
     [766] "363.1657_3.6926"                       
     [767] "363.1658_2.8504"                       
     [768] "363.1659_2.2639"                       
     [769] "363.1659_2.4664"                       
     [770] "364.14_3.0751"                         
     [771] "365.0525_0.4907"                       
     [772] "365.0643_2.2828"                       
     [773] "365.0644_1.828"                        
     [774] "365.0717_0.326"                        
     [775] "365.1087_1.7566"                       
     [776] "365.1452_3.0016"                       
     [777] "365.1453_2.1547"                       
     [778] "365.1579_2.1532"                       
     [779] "365.158_2.2607"                        
     [780] "365.1582_3.0027"                       
     [781] "365.1815_2.6256"                       
     [782] "365.1815_3.585"                        
     [783] "365.1969_5.5496"                       
     [784] "366.0829_1.814"                        
     [785] "366.1403_1.1009"                       
     [786] "367.049_2.4251"                        
     [787] "367.0621_2.5156"                       
     [788] "367.0668_2.0046"                       
     [789] "367.0873_1.5748"                       
     [790] "367.0879_1.1059"                       
     [791] "367.1033_2.5713"                       
     [792] "Feruloylquinic acid"                   
     [793] "367.1163_2.2899"                       
     [794] "367.1243_1.7697"                       
     [795] "367.1606_2.1812"                       
     [796] "367.1606_2.8802"                       
     [797] "367.2487_5.8805"                       
     [798] "368.0988_1.7827"                       
     [799] "368.1112_1.7833"                       
     [800] "369.0066_0.4796"                       
     [801] "369.0198_0.4777"                       
     [802] "369.0825_2.4095"                       
     [803] "369.0835_2.0896"                       
     [804] "369.1031_0.5005"                       
     [805] "369.1342_5.0371"                       
     [806] "369.1764_1.9013"                       
     [807] "369.1764_2.7261"                       
     [808] "370.1039_1.7824"                       
     [809] "370.9689_0.8591"                       
     [810] "370.983_0.4807"                        
     [811] "371.0982_2.0161"                       
     [812] "371.0982_2.5947"                       
     [813] "371.0983_1.6923"                       
     [814] "371.1113_1.6917"                       
     [815] "371.1188_0.3093"                       
     [816] "371.1345_2.2679"                       
     [817] "371.1499_4.9296"                       
     [818] "371.1499_4.1926"                       
     [819] "373.0774_1.7922"                       
     [820] "373.0775_2.1389"                       
     [821] "373.1138_2.0659"                       
     [822] "373.1138_1.921"                        
     [823] "373.1502_3.5313"                       
     [824] "373.1656_3.8454"                       
     [825] "373.1865_3.3416"                       
     [826] "373.1866_2.9893"                       
     [827] "373.1866_2.8159"                       
     [828] "374.1566_0.9483"                       
     [829] "374.9599_4.3219"                       
     [830] "375.0697_2.2189"                       
     [831] "375.0698_1.9725"                       
     [832] "375.0931_1.7885"                       
     [833] "375.106_1.6253"                        
     [834] "375.1306_2.4225"                       
     [835] "375.1698_3.8434"                       
     [836] "377.0094_0.8162"                       
     [837] "377.0078_0.8054"                       
     [838] "377.0402_1.8047"                       
     [839] "377.0531_1.1015"                       
     [840] "377.0849_0.3459"                       
     [841] "377.0853_1.9986"                       
     [842] "377.0854_0.7621"                       
     [843] "377.0977_0.3693"                       
     [844] "377.1006_2.4188"                       
     [845] "377.1036_1.6255"                       
     [846] "377.1815_2.5497"                       
     [847] "377.1816_2.6991"                       
     [848] "377.1947_2.6989"                       
     [849] "378.0883_0.3882"                       
     [850] "378.1339_1.748"                        
     [851] "378.9597_0.4279"                       
     [852] "379.0174_0.8566"                       
     [853] "379.0192_1.1024"                       
     [854] "379.0556_1.8037"                       
     [855] "379.0817_0.3199"                       
     [856] "379.0832_1.5242"                       
     [857] "379.0983_0.3112"                       
     [858] "379.1244_1.9887"                       
     [859] "379.1608_1.914"                        
     [860] "379.1608_2.2789"                       
     [861] "379.1608_2.5857"                       
     [862] "379.1971_2.5498"                       
     [863] "379.1971_2.9338"                       
     [864] "380.0444_2.466"                        
     [865] "380.1561_1.8479"                       
     [866] "381.0615_1.8035"                       
     [867] "381.1319_2.6639"                       
     [868] "381.1402_1.942"                        
     [869] "381.1531_1.9423"                       
     [870] "381.1765_2.5099"                       
     [871] "381.1896_2.5097"                       
     [872] "381.2127_3.0704"                       
     [873] "382.9835_1.8038"                       
     [874] "383.1192_1.1066"                       
     [875] "383.1293_2.6636"                       
     [876] "384.9807_0.4791"                       
     [877] "384.9937_0.4777"                       
     [878] "385.0563_2.6124"                       
     [879] "385.0596_2.0504"                       
     [880] "385.114_2.3865"                        
     [881] "385.117_1.7377"                        
     [882] "385.1269_2.3878"                       
     [883] "385.1271_2.7069"                       
     [884] "385.1349_1.625"                        
     [885] "385.1501_3.1154"                       
     [886] "385.163_2.2347"                        
     [887] "385.2149_6.3295"                       
     [888] "386.9804_0.4769"                       
     [889] "387.0033_0.4592"                       
     [890] "387.0171_0.4507"                       
     [891] "387.0593_2.869"                        
     [892] "387.0931_2.1413"                       
     [893] "387.1143_1.4992"                       
     [894] "387.1267_0.4355"                       
     [895] "387.1293_2.7863"                       
     [896] "387.1296_2.1184"                       
     [897] "387.1425_2.1181"                       
     [898] "387.1447_3.5225"                       
     [899] "387.1448_3.3138"                       
     [900] "387.1656_2.4023"                       
     [901] "387.1657_3.023"                        
     [902] "388.9937_0.4779"                       
     [903] "389.1216_1.8576"                       
     [904] "389.145_2.4001"                        
     [905] "389.1451_2.1766"                       
     [906] "389.1603_3.4491"                       
     [907] "389.1604_3.3421"                       
     [908] "389.1814_2.4772"                       
     [909] "390.0355_0.437"                        
     [910] "390.0486_0.4518"                       
     [911] "390.0718_0.3377"                       
     [912] "390.9547_3.4"                          
     [913] "391.0321_0.4475"                       
     [914] "391.0435_1.9732"                       
     [915] "391.0435_2.2153"                       
     [916] "391.0645_1.8043"                       
     [917] "391.0691_0.5426"                       
     [918] "391.08_2.1267"                         
     [919] "391.1242_2"                            
     [920] "391.1607_2.3826"                       
     [921] "391.1607_1.9065"                       
     [922] "391.1608_2.6636"                       
     [923] "391.176_3.2728"                        
     [924] "391.1971_3.1865"                       
     [925] "391.9681_0.8253"                       
     [926] "392.0354_0.42"                         
     [927] "392.0487_0.5929"                       
     [928] "392.0489_0.3971"                       
     [929] "392.1349_3.3265"                       
     [930] "392.9758_0.8264"                       
     [931] "393.0123_1.8046"                       
     [932] "393.0368_0.4471"                       
     [933] "393.0368_0.3957"                       
     [934] "393.037_0.6054"                        
     [935] "393.0494_0.59"                         
     [936] "393.0743_0.3357"                       
     [937] "393.0774_2.5188"                       
     [938] "393.1033_1.6833"                       
     [939] "393.1402_2.241"                        
     [940] "393.1533_2.2405"                       
     [941] "393.1761_2.685"                        
     [942] "393.1763_1.8766"                       
     [943] "393.1764_2.0814"                       
     [944] "393.1764_1.7781"                       
     [945] "393.1764_2.9477"                       
     [946] "394.04_0.6052"                         
     [947] "395.0636_2.2126"                       
     [948] "395.1556_2.0539"                       
     [949] "395.1557_2.2"                          
     [950] "395.192_1.9342"                        
     [951] "395.192_2.9414"                        
     [952] "395.1921_2.5833"                       
     [953] "395.1921_2.2289"                       
     [954] "396.1955_2.9337"                       
     [955] "397.2077_2.3909"                       
     [956] "399.0992_2.1189"                       
     [957] "399.1122_0.3045"                       
     [958] "399.1506_1.8577"                       
     [959] "399.9362_0.4792"                       
     [960] "400.944_0.4805"                        
     [961] "400.9576_0.479"                        
     [962] "401.1087_2.1635"                       
     [963] "401.1451_2.3659"                       
     [964] "403.0046_0.4439"                       
     [965] "403.1244_1.9329"                       
     [966] "403.1374_2.1141"                       
     [967] "403.1606_2.1372"                       
     [968] "403.176_3.7783"                        
     [969] "403.197_3.8146"                        
     [970] "403.1971_4.0731"                       
     [971] "404.9911_0.4463"                       
     [972] "405.0282_1.5665"                       
     [973] "405.08_2.4616"                         
     [974] "405.14_1.6423"                         
     [975] "405.1549_2.9222"                       
     [976] "405.1765_2.8865"                       
     [977] "405.1899_2.8864"                       
     [978] "405.1917_3.5604"                       
     [979] "406.9614_0.4715"                       
     [980] "407.0069_2.2126"                       
     [981] "407.0167_0.6497"                       
     [982] "407.1557_1.9261"                       
     [983] "407.1557_2.7363"                       
     [984] "407.1557_2.4628"                       
     [985] "407.1692_2.4628"                       
     [986] "407.192_2.9992"                        
     [987] "407.192_3.4213"                        
     [988] "408.0878_2.2878"                       
     [989] "408.9449_0.4789"                       
     [990] "409.0445_1.7375"                       
     [991] "409.0568_1.8026"                       
     [992] "409.1502_3.2205"                       
     [993] "409.1712_2.1498"                       
     [994] "409.2076_2.9911"                       
     [995] "411.0461_2.2199"                       
     [996] "411.0917_0.3098"                       
     [997] "411.1869_3.0404"                       
     [998] "411.1869_2.2889"                       
     [999] "411.1869_1.9994"                       
    [1000] "413.1663_2.1023"                       
    [1001] "413.2178_4.2441"                       
    [1002] "413.2179_3.2994"                       
    [1003] "413.2444_5.7394"                       
    [1004] "414.9592_0.8284"                       
    [1005] "415.1607_2.6307"                       
    [1006] "415.1971_3.1078"                       
    [1007] "416.0833_2.2207"                       
    [1008] "417.0823_2.3829"                       
    [1009] "417.0824_3.4837"                       
    [1010] "417.0825_3.1673"                       
    [1011] "417.1037_2.1755"                       
    [1012] "417.14_2.2004"                         
    [1013] "417.1532_2.4248"                       
    [1014] "417.1668_2.4286"                       
    [1015] "417.1763_3.748"                        
    [1016] "417.2127_3.146"                        
    [1017] "417.9525_1.8041"                       
    [1018] "418.9545_0.4527"                       
    [1019] "418.9601_1.8031"                       
    [1020] "419.0439_1.5459"                       
    [1021] "419.0805_4.825"                        
    [1022] "419.1164_2.4807"                       
    [1023] "419.1506_2.4249"                       
    [1024] "419.192_3.3221"                        
    [1025] "419.192_2.816"                         
    [1026] "419.1921_2.9878"                       
    [1027] "419.2064_3.3339"                       
    [1028] "420.0561_4.6306"                       
    [1029] "420.1954_2.8159"                       
    [1030] "420.1954_2.9877"                       
    [1031] "420.2024_2.5743"                       
    [1032] "420.9951_1.569"                        
    [1033] "421.0007_0.6693"                       
    [1034] "421.0018_0.4571"                       
    [1035] "421.0155_0.4798"                       
    [1036] "421.0745_2.2172"                       
    [1037] "421.1348_1.7046"                       
    [1038] "421.135_2.3849"                        
    [1039] "421.1634_2.4554"                       
    [1040] "421.1702_2.6222"                       
    [1041] "421.1714_2.2182"                       
    [1042] "421.2075_3.1891"                       
    [1043] "421.2077_3.6193"                       
    [1044] "421.2217_3.6193"                       
    [1045] "422.005_0.4696"                        
    [1046] "422.0239_0.5153"                       
    [1047] "422.0531_4.6311"                       
    [1048] "422.1667_2.456"                        
    [1049] "422.9354_0.4706"                       
    [1050] "423.0018_0.4733"                       
    [1051] "423.0598_1.7892"                       
    [1052] "423.1606_2.456"                        
    [1053] "423.1786_2.4766"                       
    [1054] "423.1791_2.6652"                       
    [1055] "423.1862_2.4094"                       
    [1056] "423.2233_3.4312"                       
    [1057] "423.2234_3.3419"                       
    [1058] "424.9719_0.4087"                       
    [1059] "425.1759_2.4776"                       
    [1060] "425.1814_4.279"                        
    [1061] "425.1814_2.8927"                       
    [1062] "425.1966_3.791"                        
    [1063] "425.2026_2.7214"                       
    [1064] "426.9544_1.8036"                       
    [1065] "427.1455_1.9423"                       
    [1066] "427.1607_3.0629"                       
    [1067] "427.1608_2.6028"                       
    [1068] "427.1824_2.424"                        
    [1069] "427.1956_2.4889"                       
    [1070] "427.1957_2.2824"                       
    [1071] "427.1969_2.6918"                       
    [1072] "428.1991_2.4271"                       
    [1073] "428.9515_1.8038"                       
    [1074] "429.1247_1.1065"                       
    [1075] "429.153_2.3659"                        
    [1076] "429.1611_1.6565"                       
    [1077] "429.1764_3.2774"                       
    [1078] "429.1765_3.0628"                       
    [1079] "429.1901_3.0642"                       
    [1080] "429.1906_3.277"                        
    [1081] "429.2127_2.9865"                       
    [1082] "429.2128_3.3441"                       
    [1083] "430.9886_0.8559"                       
    [1084] "431.0651_2.1848"                       
    [1085] "431.0982_2.7665"                       
    [1086] "431.1123_2.7677"                       
    [1087] "431.1397_0.4992"                       
    [1088] "431.1688_2.5839"                       
    [1089] "431.1688_2.8646"                       
    [1090] "431.1708_4.3389"                       
    [1091] "431.1709_4.1131"                       
    [1092] "431.1827_2.8668"                       
    [1093] "431.1913_3.12"                         
    [1094] "431.192_3.3155"                        
    [1095] "431.1922_2.4561"                       
    [1096] "431.2058_2.4541"                       
    [1097] "431.2203_5.9789"                       
    [1098] "431.2283_3.3653"                       
    [1099] "433.0411_2.6812"                       
    [1100] "433.0551_2.6794"                       
    [1101] "433.079_2.1223"                        
    [1102] "433.0911_2.9514"                       
    [1103] "433.1138_3.0804"                       
    [1104] "433.1284_2.1683"                       
    [1105] "433.135_1.5891"                        
    [1106] "433.1477_1.7816"                       
    [1107] "433.1663_2.5834"                       
    [1108] "433.1667_2.8641"                       
    [1109] "433.1865_3.9498"                       
    [1110] "433.2075_2.4761"                       
    [1111] "433.2075_3.3919"                       
    [1112] "433.2076_2.6656"                       
    [1113] "433.2077_2.7858"                       
    [1114] "433.2216_2.4767"                       
    [1115] "433.2229_3.7635"                       
    [1116] "433.2357_6.4562"                       
    [1117] "434.0949_2.9672"                       
    [1118] "434.211_2.4755"                        
    [1119] "434.211_2.6646"                        
    [1120] "435.0109_1.5463"                       
    [1121] "435.0359_1.671"                        
    [1122] "435.083_2.9706"                        
    [1123] "435.0838_2.4801"                       
    [1124] "435.0929_2.1316"                       
    [1125] "435.1081_2.0985"                       
    [1126] "Phlorhizin"                            
    [1127] "435.1505_1.7206"                       
    [1128] "435.166_2.2206"                        
    [1129] "435.1869_2.9512"                       
    [1130] "435.223_2.547"                         
    [1131] "435.2232_3.4991"                       
    [1132] "436.1468_3.1803"                       
    [1133] "437.0173_2.463"                        
    [1134] "437.04_2.3985"                         
    [1135] "437.0401_2.2184"                       
    [1136] "437.0958_0.3118"                       
    [1137] "437.1086_2.0042"                       
    [1138] "437.1087_2.5655"                       
    [1139] "437.1089_2.2905"                       
    [1140] "437.1299_2.2021"                       
    [1141] "437.1662_2.3646"                       
    [1142] "438.0435_2.2163"                       
    [1143] "438.8998_0.4724"                       
    [1144] "439.0034_0.404"                        
    [1145] "439.055_1.8238"                        
    [1146] "439.0829_0.3097"                       
    [1147] "439.1229_0.3102"                       
    [1148] "439.1814_2.6228"                       
    [1149] "439.2182_2.831"                        
    [1150] "440.0841_0.3217"                       
    [1151] "440.0872_0.3084"                       
    [1152] "440.9467_0.3803"                       
    [1153] "441.0255_2.4804"                       
    [1154] "441.0901_0.3136"                       
    [1155] "441.1222_2.8658"                       
    [1156] "441.1764_2.7634"                       
    [1157] "441.1979_2.5839"                       
    [1158] "441.198_2.8652"                        
    [1159] "441.2114_2.9115"                       
    [1160] "441.2415_2.8652"                       
    [1161] "442.215_2.8693"                        
    [1162] "442.9822_0.4322"                       
    [1163] "442.9835_2.2251"                       
    [1164] "442.9835_2.3985"                       
    [1165] "443.1555_3.4521"                       
    [1166] "443.1768_1.9931"                       
    [1167] "443.1768_1.7812"                       
    [1168] "443.192_3.0964"                        
    [1169] "443.1921_2.1252"                       
    [1170] "443.2283_4.5586"                       
    [1171] "444.1719_2.439"                        
    [1172] "444.9649_0.5037"                       
    [1173] "444.9665_0.4387"                       
    [1174] "445.1351_1.8287"                       
    [1175] "445.1713_2.3354"                       
    [1176] "445.1843_2.9911"                       
    [1177] "445.1863_4.5421"                       
    [1178] "445.2076_3.3403"                       
    [1179] "446.9628_0.6403"                       
    [1180] "446.9644_0.4719"                       
    [1181] "447.0718_3.5779"                       
    [1182] "447.093_2.5675"                        
    [1183] "Quercitrin"                            
    [1184] "447.107_2.1741"                        
    [1185] "447.1072_3.0181"                       
    [1186] "447.1142_1.9486"                       
    [1187] "447.1616_2.289"                        
    [1188] "447.1659_3.326"                        
    [1189] "447.1819_2.9915"                       
    [1190] "447.1867_2.5444"                       
    [1191] "447.2233_3.6381"                       
    [1192] "447.2376_2.8847"                       
    [1193] "448.0101_2.2137"                       
    [1194] "448.0964_2.5688"                       
    [1195] "448.0965_2.1825"                       
    [1196] "448.1108_3.0185"                       
    [1197] "448.1823_2.4561"                       
    [1198] "449.0258_2.4785"                       
    [1199] "449.1081_2.4588"                       
    [1200] "449.1087_3.1113"                       
    [1201] "449.2025_2.8009"                       
    [1202] "449.9861_6.1719"                       
    [1203] "449.9862_5.049"                        
    [1204] "450.9813_1.5457"                       
    [1205] "450.9859_2.2154"                       
    [1206] "451.0007_2.2141"                       
    [1207] "451.0403_0.8117"                       
    [1208] "451.0542_2.4796"                       
    [1209] "451.0878_2.314"                        
    [1210] "451.1032_3.2906"                       
    [1211] "451.1032_2.8192"                       
    [1212] "451.1244_2.0438"                       
    [1213] "451.1387_2.0455"                       
    [1214] "451.1407_2.4505"                       
    [1215] "451.2182_2.1765"                       
    [1216] "451.9831_6.1722"                       
    [1217] "451.9894_2.2152"                       
    [1218] "452.0213_4.7886"                       
    [1219] "452.0516_1.7813"                       
    [1220] "452.278_5.9026"                        
    [1221] "452.9838_2.2153"                       
    [1222] "453.012_2.2186"                        
    [1223] "453.0123_1.973"                        
    [1224] "453.1763_3.1049"                       
    [1225] "453.1975_2.7422"                       
    [1226] "453.9869_2.213"                        
    [1227] "454.0173_4.7928"                       
    [1228] "454.9813_2.2139"                       
    [1229] "455.1556_2.7208"                       
    [1230] "455.213_3.2738"                        
    [1231] "455.2134_2.9913"                       
    [1232] "455.3526_7.3681"                       
    [1233] "456.356_7.3686"                        
    [1234] "456.9104_0.4456"                       
    [1235] "457.0216_2.4617"                       
    [1236] "457.03_1.7886"                         
    [1237] "457.117_2.8981"                        
    [1238] "457.1192_0.3268"                       
    [1239] "457.1196_1.0996"                       
    [1240] "457.1349_2.2898"                       
    [1241] "457.1713_2.1528"                       
    [1242] "457.1865_3.6873"                       
    [1243] "457.1923_2.048"                        
    [1244] "457.1924_2.2889"                       
    [1245] "457.2076_3.3767"                       
    [1246] "458.0298_4.6314"                       
    [1247] "458.1504_0.4856"                       
    [1248] "458.1875_2.8643"                       
    [1249] "458.9509_0.6992"                       
    [1250] "458.9566_0.5285"                       
    [1251] "458.9573_0.4294"                       
    [1252] "459.0448_1.8289"                       
    [1253] "459.1345_0.7782"                       
    [1254] "459.1504_2.467"                        
    [1255] "459.1717_1.584"                        
    [1256] "459.2232_3.0714"                       
    [1257] "461.1661_2.2368"                       
    [1258] "461.1666_1.7995"                       
    [1259] "461.179_2.1156"                        
    [1260] "461.1796_2.6853"                       
    [1261] "461.1806_2.4556"                       
    [1262] "461.1813_3.174"                        
    [1263] "461.1814_4.135"                        
    [1264] "461.2387_3.375"                        
    [1265] "462.1268_2.883"                        
    [1266] "462.9652_2.4613"                       
    [1267] "Quercetin-3-glucoside"                 
    [1268] "463.0884_2.3069"                       
    [1269] "463.0886_2.4747"                       
    [1270] "463.0886_2.012"                        
    [1271] "463.1023_2.8037"                       
    [1272] "463.1235_2.9076"                       
    [1273] "463.1559_2.8921"                       
    [1274] "463.1818_2.2781"                       
    [1275] "463.1968_2.6056"                       
    [1276] "463.2177_2.8115"                       
    [1277] "464.1035_2.3208"                       
    [1278] "464.1061_2.7947"                       
    [1279] "464.9623_2.4613"                       
    [1280] "464.9969_1.5402"                       
    [1281] "465.0732_0.856"                        
    [1282] "465.1037_2.6907"                       
    [1283] "465.1527_2.6915"                       
    [1284] "465.1528_3.0629"                       
    [1285] "465.1968_3.1644"                       
    [1286] "466.1213_2.1741"                       
    [1287] "466.1214_2.2892"                       
    [1288] "467.0352_0.8425"                       
    [1289] "467.0361_2.5157"                       
    [1290] "467.0717_2.1283"                       
    [1291] "467.0805_2.2206"                       
    [1292] "467.1164_0.3286"                       
    [1293] "467.1192_1.7544"                       
    [1294] "467.1386_1.9608"                       
    [1295] "467.1708_4.572"                        
    [1296] "468.1508_2.3554"                       
    [1297] "469.0806_2.6396"                       
    [1298] "469.0905_3.597"                        
    [1299] "469.171_2.4292"                        
    [1300] "469.1713_3.5302"                       
    [1301] "469.1865_3.7907"                       
    [1302] "469.2011_4.5771"                       
    [1303] "469.2287_3.6169"                       
    [1304] "469.3319_7.0007"                       
    [1305] "470.3353_6.9983"                       
    [1306] "471.1345_0.3133"                       
    [1307] "471.1716_1.6165"                       
    [1308] "471.1869_2.8908"                       
    [1309] "471.208_2.115"                         
    [1310] "471.2081_2.7207"                       
    [1311] "471.2235_3.4108"                       
    [1312] "473.0021_1.79"                         
    [1313] "473.1506_0.5041"                       
    [1314] "473.1509_0.3245"                       
    [1315] "473.1659_3.1674"                       
    [1316] "473.166_2.8134"                        
    [1317] "473.1814_3.8139"                       
    [1318] "474.1195_2.0108"                       
    [1319] "474.1693_3.1674"                       
    [1320] "474.2623_5.2767"                       
    [1321] "474.9209_0.4396"                       
    [1322] "475.1455_2.1806"                       
    [1323] "475.1456_1.9221"                       
    [1324] "475.1487_1.6154"                       
    [1325] "475.1605_2.8874"                       
    [1326] "475.1949_2.8271"                       
    [1327] "475.1971_4.5212"                       
    [1328] "475.2182_2.9064"                       
    [1329] "475.2182_3.3694"                       
    [1330] "476.2215_2.9061"                       
    [1331] "476.278_5.59"                          
    [1332] "477.0021_2.4807"                       
    [1333] "477.1035_3.0035"                       
    [1334] "477.1154_1.6445"                       
    [1335] "477.1246_1.7546"                       
    [1336] "477.161_2.1769"                        
    [1337] "477.1927_2.8255"                       
    [1338] "477.2111_2.8327"                       
    [1339] "477.233_2.8549"                        
    [1340] "477.2493_4.4304"                       
    [1341] "478.099_1.7162"                        
    [1342] "478.2935_5.9807"                       
    [1343] "478.9601_2.233"                        
    [1344] "478.9992_2.4804"                       
    [1345] "479.1192_2.6422"                       
    [1346] "479.1531_2.4237"                       
    [1347] "479.1768_3.3084"                       
    [1348] "480.3092_5.9362"                       
    [1349] "481.0509_0.857"                        
    [1350] "481.0655_0.8528"                       
    [1351] "481.0879_2.4805"                       
    [1352] "481.0982_2.2996"                       
    [1353] "481.1484_2.3342"                       
    [1354] "481.1843_3.2547"                       
    [1355] "481.1843_2.8276"                       
    [1356] "482.0484_5.585"                        
    [1357] "483.0454_2.2157"                       
    [1358] "483.0697_2.6201"                       
    [1359] "483.1716_2.0146"                       
    [1360] "483.2002_3.4811"                       
    [1361] "483.2081_2.5438"                       
    [1362] "483.2724_6.681"                        
    [1363] "483.2724_6.9429"                       
    [1364] "483.3113_6.1796"                       
    [1365] "484.2677_4.7064"                       
    [1366] "484.9963_2.4787"                       
    [1367] "485.0611_2.3671"                       
    [1368] "485.0676_3.0162"                       
    [1369] "485.1507_1.6193"                       
    [1370] "485.1661_2.9328"                       
    [1371] "485.1662_2.5973"                       
    [1372] "485.1775_2.3517"                       
    [1373] "485.2025_2.4383"                       
    [1374] "485.2231_2.4834"                       
    [1375] "485.2238_2.8252"                       
    [1376] "485.2383_2.8185"                       
    [1377] "485.2389_4.903"                        
    [1378] "485.3269_4.7447"                       
    [1379] "485.3271_5.8609"                       
    [1380] "486.0427_5.5835"                       
    [1381] "486.3454_5.8658"                       
    [1382] "486.9545_2.2306"                       
    [1383] "486.9933_2.4787"                       
    [1384] "487.101_2.9602"                        
    [1385] "487.1958_2.8099"                       
    [1386] "487.197_4.167"                         
    [1387] "487.3426_5.579"                        
    [1388] "487.3426_5.0473"                       
    [1389] "487.3427_5.2334"                       
    [1390] "487.3574_5.237"                        
    [1391] "487.3575_5.0493"                       
    [1392] "488.9515_2.2299"                       
    [1393] "488.9912_2.4789"                       
    [1394] "489.1247_2.0714"                       
    [1395] "489.1457_1.1002"                       
    [1396] "489.1975_3.1759"                       
    [1397] "489.2126_3.8329"                       
    [1398] "490.1018_2.9602"                       
    [1399] "490.9494_2.2294"                       
    [1400] "490.9893_2.4791"                       
    [1401] "491.0464_2.7117"                       
    [1402] "491.1767_2.3942"                       
    [1403] "491.1902_2.6298"                       
    [1404] "491.1918_3.7398"                       
    [1405] "491.1919_4.0801"                       
    [1406] "491.1921_3.8131"                       
    [1407] "491.2126_2.827"                        
    [1408] "491.213_2.44"                          
    [1409] "491.2271_3.114"                        
    [1410] "491.2282_3.425"                        
    [1411] "492.0771_5.5841"                       
    [1412] "492.1276_2.3036"                       
    [1413] "492.2159_2.8239"                       
    [1414] "492.2164_3.252"                        
    [1415] "493.1195_1.6163"                       
    [1416] "493.1348_2.9036"                       
    [1417] "493.1687_2.8651"                       
    [1418] "493.1688_2.5841"                       
    [1419] "493.2286_3.6469"                       
    [1420] "493.2287_3.4805"                       
    [1421] "493.244_3.4855"                        
    [1422] "494.1251_2.3032"                       
    [1423] "494.2321_3.4806"                       
    [1424] "495.1108_0.3078"                       
    [1425] "495.1163_2.8653"                       
    [1426] "495.2444_3.9488"                       
    [1427] "497.009_0.4328"                        
    [1428] "497.0231_0.8565"                       
    [1429] "497.1222_2.0098"                       
    [1430] "497.1871_2.0563"                       
    [1431] "497.1872_2.3587"                       
    [1432] "497.1872_2.2513"                       
    [1433] "497.2601_4.0878"                       
    [1434] "499.0117_2.1564"                       
    [1435] "499.1946_2.8631"                       
    [1436] "499.2546_3.949"                        
    [1437] "501.0096_2.1565"                       
    [1438] "501.1307_2.0958"                       
    [1439] "501.1457_0.7657"                       
    [1440] "501.1457_1.1056"                       
    [1441] "501.1458_1.5295"                       
    [1442] "501.1763_4.0936"                       
    [1443] "501.1763_3.6653"                       
    [1444] "501.2185_2.6283"                       
    [1445] "501.2702_5.3958"                       
    [1446] "501.3217_3.6013"                       
    [1447] "501.3218_4.5814"                       
    [1448] "501.3218_4.865"                        
    [1449] "501.3219_5.11"                         
    [1450] "501.3219_5.4605"                       
    [1451] "502.1565_2.3029"                       
    [1452] "503.1402_2.0455"                       
    [1453] "503.1611_0.504"                        
    [1454] "503.2027_3.2122"                       
    [1455] "503.2494_3.6076"                       
    [1456] "503.2495_3.9548"                       
    [1457] "503.3368_4.8437"                       
    [1458] "503.3374_3.7115"                       
    [1459] "503.3379_4.2797"                       
    [1460] "504.1876_3.2298"                       
    [1461] "504.3092_5.6198"                       
    [1462] "504.3406_4.8429"                       
    [1463] "505.0934_3.1768"                       
    [1464] "505.0985_2.9161"                       
    [1465] "505.1196_1.6183"                       
    [1466] "505.1564_1.7774"                       
    [1467] "505.2568_5.8111"                       
    [1468] "505.2592_6.1931"                       
    [1469] "505.3593_4.2848"                       
    [1470] "506.1512_1.5499"                       
    [1471] "507.0375_0.4762"                       
    [1472] "507.1141_2.4349"                       
    [1473] "507.1352_1.6321"                       
    [1474] "507.1715_1.7987"                       
    [1475] "507.2079_3.6713"                       
    [1476] "507.2232_3.128"                        
    [1477] "507.2724_6.3168"                       
    [1478] "507.9263_0.6126"                       
    [1479] "508.9255_0.6109"                       
    [1480] "509.1508_1.8286"                       
    [1481] "509.1872_2.5973"                       
    [1482] "509.2227_2.8685"                       
    [1483] "510.0796_6.329"                        
    [1484] "510.226_2.8694"                        
    [1485] "510.926_0.6077"                        
    [1486] "511.1817_3.7069"                       
    [1487] "511.1818_2.9995"                       
    [1488] "512.0768_6.3285"                       
    [1489] "512.299_5.2989"                        
    [1490] "513.1111_2.8053"                       
    [1491] "513.1458_1.5163"                       
    [1492] "513.1458_1.1047"                       
    [1493] "513.2545_2.9176"                       
    [1494] "514.074_6.3285"                        
    [1495] "514.1827_2.3922"                       
    [1496] "515.1194_2.9832"                       
    [1497] "515.1247_0.3234"                       
    [1498] "515.1251_1.1032"                       
    [1499] "515.1251_1.5634"                       
    [1500] "515.1402_2.0627"                       
    [1501] "515.1402_1.7862"                       
    [1502] "515.1405_2.4686"                       
    [1503] "515.1615_1.3251"                       
    [1504] "515.1615_1.5651"                       
    [1505] "515.1767_2.4707"                       
    [1506] "515.1919_3.8105"                       
    [1507] "515.1927_4.5717"                       
    [1508] "515.2124_3.4635"                       
    [1509] "515.2342_2.7016"                       
    [1510] "516.2377_2.7016"                       
    [1511] "517.1194_1.7956"                       
    [1512] "517.1399_0.3081"                       
    [1513] "517.1557_2.4914"                       
    [1514] "517.1557_2.0082"                       
    [1515] "517.1587_0.309"                        
    [1516] "517.1707_3.291"                        
    [1517] "517.2282_3.3109"                       
    [1518] "517.2286_3.0066"                       
    [1519] "517.2287_3.5616"                       
    [1520] "517.2288_2.4344"                       
    [1521] "517.3167_3.9765"                       
    [1522] "517.3168_4.2139"                       
    [1523] "518.3199_3.9765"                       
    [1524] "518.9655_0.533"                        
    [1525] "518.9931_0.4686"                       
    [1526] "519.003_0.855"                         
    [1527] "519.1465_2.3042"                       
    [1528] "519.1867_3.3225"                       
    [1529] "519.1868_2.4973"                       
    [1530] "519.2231_3.8992"                       
    [1531] "519.2236_3.3612"                       
    [1532] "519.2443_2.5828"                       
    [1533] "519.2443_3.4473"                       
    [1534] "519.3323_3.5282"                       
    [1535] "519.3323_4.0485"                       
    [1536] "519.3324_4.5605"                       
    [1537] "520.2677_5.6714"                       
    [1538] "521.2014_2.5235"                       
    [1539] "521.2024_2.3254"                       
    [1540] "521.2026_2.7914"                       
    [1541] "521.2389_3.2981"                       
    [1542] "521.2391_3.3992"                       
    [1543] "522.1092_1.9851"                       
    [1544] "522.1824_0.9113"                       
    [1545] "522.9756_2.3871"                       
    [1546] "523.1035_3.177"                        
    [1547] "523.203_2.5634"                        
    [1548] "523.2177_2.6659"                       
    [1549] "523.2179_2.9759"                       
    [1550] "523.34_7.3716"                         
    [1551] "525.0251_3.1773"                       
    [1552] "525.0407_0.8449"                       
    [1553] "525.1821_2.2991"                       
    [1554] "525.2182_2.5811"                       
    [1555] "526.0441_0.8396"                       
    [1556] "527.104_1.7488"                        
    [1557] "527.1402_2.532"                        
    [1558] "528.0412_0.4037"                       
    [1559] "529.156_2.2033"                        
    [1560] "529.2288_2.7738"                       
    [1561] "530.3015_5.935"                        
    [1562] "531.1503_3.1126"                       
    [1563] "531.1867_3.7944"                       
    [1564] "531.1868_3.5627"                       
    [1565] "531.2229_3.3903"                       
    [1566] "531.2231_3.6509"                       
    [1567] "531.2233_4.6026"                       
    [1568] "531.2807_5.1519"                       
    [1569] "532.1821_4.572"                        
    [1570] "533.1347_0.3092"                       
    [1571] "533.1656_3.0548"                       
    [1572] "533.1716_0.3161"                       
    [1573] "533.2018_3.6771"                       
    [1574] "533.2022_2.8782"                       
    [1575] "533.2234_2.6546"                       
    [1576] "533.3116_4.135"                        
    [1577] "534.9707_0.857"                        
    [1578] "535.0245_3.1766"                       
    [1579] "535.0538_3.177"                        
    [1580] "535.122_2.0909"                        
    [1581] "535.2391_3.4237"                       
    [1582] "537.1056_2.5169"                       
    [1583] "537.1246_3.1368"                       
    [1584] "537.1457_1.705"                        
    [1585] "537.1616_2.3741"                       
    [1586] "537.1651_0.3088"                       
    [1587] "537.1974_3.2839"                       
    [1588] "538.3146_5.4955"                       
    [1589] "539.1872_4.6053"                       
    [1590] "539.1977_2.521"                        
    [1591] "539.2342_2.2309"                       
    [1592] "539.3141_4.4759"                       
    [1593] "539.3298_4.2868"                       
    [1594] "540.3305_5.9353"                       
    [1595] "541.1843_3.0854"                       
    [1596] "541.3121_4.2793"                       
    [1597] "542.1776_2.7418"                       
    [1598] "543.0106_0.4807"                       
    [1599] "543.1845_2.4633"                       
    [1600] "543.1999_3.1281"                       
    [1601] "545.1793_2.8949"                       
    [1602] "545.2025_3.9912"                       
    [1603] "547.0406_2.7923"                       
    [1604] "547.1665_1.9672"                       
    [1605] "547.1665_2.2008"                       
    [1606] "547.1773_2.8947"                       
    [1607] "547.2179_3.8807"                       
    [1608] "547.2181_4.606"                        
    [1609] "547.2182_3.4662"                       
    [1610] "547.2392_2.3161"                       
    [1611] "549.0834_2.2153"                       
    [1612] "549.1083_2.2113"                       
    [1613] "549.1457_1.7783"                       
    [1614] "549.1893_4.8377"                       
    [1615] "549.1951_2.3031"                       
    [1616] "549.1976_3.2211"                       
    [1617] "549.2549_2.173"                        
    [1618] "549.6102_2.2119"                       
    [1619] "551.104_2.3864"                        
    [1620] "551.1809_0.3169"                       
    [1621] "551.213_3.0839"                        
    [1622] "551.2339_2.6241"                       
    [1623] "553.0805_0.364"                        
    [1624] "553.1005_2.3081"                       
    [1625] "553.1188_2.0478"                       
    [1626] "553.1842_3.599"                        
    [1627] "553.2057_2.4339"                       
    [1628] "553.3141_5.8609"                       
    [1629] "554.145_2.3148"                        
    [1630] "554.3014_5.7252"                       
    [1631] "554.6475_2.3146"                       
    [1632] "555.0585_3.1765"                       
    [1633] "555.1353_2.112"                        
    [1634] "555.1634_2.7747"                       
    [1635] "556.2992_5.6206"                       
    [1636] "556.317_6.1257"                        
    [1637] "557.1297_2.3442"                       
    [1638] "557.2364_2.5152"                       
    [1639] "557.2441_2.5601"                       
    [1640] "558.9882_0.4808"                       
    [1641] "559.2241_2.3024"                       
    [1642] "559.2757_4.5868"                       
    [1643] "559.3119_5.8278"                       
    [1644] "559.9943_3.1804"                       
    [1645] "560.9931_0.4806"                       
    [1646] "561.0017_3.177"                        
    [1647] "561.0876_2.0982"                       
    [1648] "561.1429_2.8112"                       
    [1649] "561.1973_3.5278"                       
    [1650] "562.1433_2.2396"                       
    [1651] "562.3146_5.4026"                       
    [1652] "562.9906_0.8414"                       
    [1653] "562.9991_3.1779"                       
    [1654] "563.1403_2.6905"                       
    [1655] "563.1613_1.6351"                       
    [1656] "563.1613_1.892"                        
    [1657] "563.213_3.6229"                        
    [1658] "563.2131_2.9751"                       
    [1659] "563.2345_2.4353"                       
    [1660] "563.2827_2.4352"                       
    [1661] "563.318_5.3126"                        
    [1662] "563.318_5.4024"                        
    [1663] "564.1436_2.6921"                       
    [1664] "564.3304_5.6196"                       
    [1665] "564.3304_5.7253"                       
    [1666] "565.0426_0.4652"                       
    [1667] "565.0464_0.3348"                       
    [1668] "565.1924_2.7741"                       
    [1669] "565.2135_2.4099"                       
    [1670] "565.2287_3.3453"                       
    [1671] "565.2289_2.9194"                       
    [1672] "565.2445_2.4356"                       
    [1673] "565.2498_2.2644"                       
    [1674] "566.3459_6.1257"                       
    [1675] "567.1723_2.9745"                       
    [1676] "567.2653_2.5128"                       
    [1677] "568.2688_2.5146"                       
    [1678] "569.002_3.1768"                        
    [1679] "569.1406_2.6919"                       
    [1680] "569.172_2.1785"                        
    [1681] "570.9994_3.1764"                       
    [1682] "571.0911_0.4087"                       
    [1683] "571.2392_2.9522"                       
    [1684] "571.2595_2.9291"                       
    [1685] "571.2886_6.277"                        
    [1686] "571.3248_4.2797"                       
    [1687] "572.9957_3.1761"                       
    [1688] "573.1035_1.6552"                       
    [1689] "573.1035_3.0854"                       
    [1690] "573.1036_2.0387"                       
    [1691] "574.2055_2.5565"                       
    [1692] "575.0141_0.3455"                       
    [1693] "575.1189_2.0025"                       
    [1694] "575.1191_2.9918"                       
    [1695] "575.2705_2.6423"                       
    [1696] "576.1258_2.8626"                       
    [1697] "576.1271_2.0993"                       
    [1698] "576.1271_2.3403"                       
    [1699] "576.1439_2.3404"                       
    [1700] "576.6457_2.5177"                       
    [1701] "576.9989_0.4572"                       
    [1702] "577.0509_3.0166"                       
    [1703] "577.1349_2.7674"                       
    [1704] "577.1352_2.0705"                       
    [1705] "Procyanidin B2"                        
    [1706] "577.1511_2.8027"                       
    [1707] "577.1564_2.9419"                       
    [1708] "577.2861_3.7247"                       
    [1709] "577.2861_4.3463"                       
    [1710] "577.6317_2.34"                         
    [1711] "577.6501_2.1929"                       
    [1712] "578.0021_0.4582"                       
    [1713] "578.1529_2.2814"                       
    [1714] "578.9599_0.8455"                       
    [1715] "579.1485_2.4261"                       
    [1716] "Naringin"                              
    [1717] "579.1759_0.3106"                       
    [1718] "579.2049_2.2872"                       
    [1719] "579.208_3.0072"                        
    [1720] "579.2654_2.6985"                       
    [1721] "580.1798_0.31"                         
    [1722] "580.2246_2.4367"                       
    [1723] "580.9734_0.4762"                       
    [1724] "581.0896_2.2156"                       
    [1725] "581.1273_2.1221"                       
    [1726] "581.1284_3.1777"                       
    [1727] "581.1507_2.3982"                       
    [1728] "581.1793_3.9491"                       
    [1729] "581.208_2.4976"                        
    [1730] "583.1431_1.9656"                       
    [1731] "583.1655_2.3962"                       
    [1732] "583.1667_2.7752"                       
    [1733] "583.2158_2.3171"                       
    [1734] "584.2059_2.3925"                       
    [1735] "584.9992_2.9597"                       
    [1736] "586.0468_0.415"                        
    [1737] "586.9969_2.9599"                       
    [1738] "587.1404_2.2788"                       
    [1739] "587.1585_0.3112"                       
    [1740] "587.2341_2.5367"                       
    [1741] "587.2988_4.2822"                       
    [1742] "589.1556_2.3331"                       
    [1743] "589.1558_3.514"                        
    [1744] "589.2344_2.2865"                       
    [1745] "591.0663_2.8833"                       
    [1746] "591.114_2.4583"                        
    [1747] "591.1141_2.1777"                       
    [1748] "591.1503_3.2229"                       
    [1749] "591.1513_2.8416"                       
    [1750] "591.1515_2.5072"                       
    [1751] "591.1562_1.7636"                       
    [1752] "591.2291_2.6762"                       
    [1753] "591.2653_2.3512"                       
    [1754] "591.2654_2.589"                        
    [1755] "591.3533_5.5576"                       
    [1756] "593.1297_2.5898"                       
    [1757] "Kampferol -3-O-glucosyl(1-2)rhamnoside"
    [1758] "593.1514_2.3876"                       
    [1759] "593.1719_1.9651"                       
    [1760] "593.1873_3.5184"                       
    [1761] "593.1925_0.5032"                       
    [1762] "593.2214_2.5599"                       
    [1763] "593.2236_3.4679"                       
    [1764] "593.2448_2.3163"                       
    [1765] "593.2728_5.4012"                       
    [1766] "594.1396_1.8955"                       
    [1767] "595.0094_0.445"                        
    [1768] "595.0098_0.6895"                       
    [1769] "595.156_2.2208"                        
    [1770] "Eriocitrin"                            
    [1771] "595.2027_2.7884"                       
    [1772] "595.2211_2.7228"                       
    [1773] "595.2604_2.173"                        
    [1774] "595.2885_5.8109"                       
    [1775] "595.2966_3.412"                        
    [1776] "595.971_3.1806"                        
    [1777] "596.0129_0.4485"                       
    [1778] "596.2201_2.7225"                       
    [1779] "596.9479_0.4748"                       
    [1780] "597.1454_2.2622"                       
    [1781] "597.1818_2.874"                        
    [1782] "597.304_6.3866"                        
    [1783] "597.9681_3.1819"                       
    [1784] "598.1751_2.7223"                       
    [1785] "598.1854_2.8725"                       
    [1786] "599.1614_1.6963"                       
    [1787] "601.1319_0.5075"                       
    [1788] "601.1362_0.312"                        
    [1789] "601.1592_2.6226"                       
    [1790] "601.1845_0.3203"                       
    [1791] "601.2861_3.9502"                       
    [1792] "602.1353_0.5017"                       
    [1793] "603.1407_0.3205"                       
    [1794] "603.1475_0.3092"                       
    [1795] "603.2503_2.5595"                       
    [1796] "604.1725_0.3151"                       
    [1797] "604.2159_2.576"                        
    [1798] "605.1465_2.9747"                       
    [1799] "605.1893_0.3277"                       
    [1800] "605.1928_0.6118"                       
    [1801] "605.2399_7.1584"                       
    [1802] "606.1877_0.3166"                       
    [1803] "607.0637_1.871"                        
    [1804] "607.145_2.3002"                        
    [1805] "607.167_3.0343"                        
    [1806] "607.2085_1.5879"                       
    [1807] "607.2169_3.0344"                       
    [1808] "607.2239_2.7642"                       
    [1809] "607.2364_2.9303"                       
    [1810] "608.1332_1.6679"                       
    [1811] "608.3176_5.9349"                       
    [1812] "609.1453_2.0481"                       
    [1813] "609.1456_2.2643"                       
    [1814] "Rutin"                                 
    [1815] "609.1847_3.0614"                       
    [1816] "609.237_3.0626"                        
    [1817] "610.236_3.0616"                        
    [1818] "611.0959_2.3547"                       
    [1819] "611.1443_1.1023"                       
    [1820] "611.1608_2.5397"                       
    [1821] "611.1949_3.2864"                       
    [1822] "611.2382_3.0612"                       
    [1823] "611.2557_2.5796"                       
    [1824] "612.1907_3.0606"                       
    [1825] "612.2606_5.935"                        
    [1826] "612.9145_0.4768"                       
    [1827] "613.1673_2.2219"                       
    [1828] "613.2134_2.6101"                       
    [1829] "613.2714_2.6681"                       
    [1830] "614.9557_0.4545"                       
    [1831] "615.2039_2.435"                        
    [1832] "615.2289_2.6439"                       
    [1833] "615.2847_5.8597"                       
    [1834] "617.1485_2.7232"                       
    [1835] "617.266_2.9296"                        
    [1836] "617.3002_5.0479"                       
    [1837] "617.3003_5.2359"                       
    [1838] "619.1246_2.7239"                       
    [1839] "620.6396_2.3717"                       
    [1840] "621.0905_2.9744"                       
    [1841] "621.1609_2.517"                        
    [1842] "623.1619_2.49"                         
    [1843] "623.1977_3.3037"                       
    [1844] "625.136_1.8046"                        
    [1845] "625.1404_2.317"                        
    [1846] "625.1407_2.6495"                       
    [1847] "626.1782_2.4755"                       
    [1848] "629.1274_2.825"                        
    [1849] "629.1283_2.9466"                       
    [1850] "629.1639_3.5192"                       
    [1851] "629.172_1.6636"                        
    [1852] "629.2083_2.8144"                       
    [1853] "631.1253_3.1758"                       
    [1854] "631.1271_2.9645"                       
    [1855] "631.1436_2.723"                        
    [1856] "631.1621_3.5181"                       
    [1857] "631.1649_3.0491"                       
    [1858] "631.2403_3.2905"                       
    [1859] "632.3174_5.6194"                       
    [1860] "632.9653_0.4388"                       
    [1861] "633.1243_2.0622"                       
    [1862] "633.1416_2.7226"                       
    [1863] "633.2365_2.5799"                       
    [1864] "633.2548_3.4058"                       
    [1865] "633.3791_6.1398"                       
    [1866] "634.9043_0.4708"                       
    [1867] "635.158_2.975"                         
    [1868] "635.172_2.2741"                        
    [1869] "635.2033_1.8047"                       
    [1870] "636.2221_2.803"                        
    [1871] "636.5316_6.7983"                       
    [1872] "637.049_3.6325"                        
    [1873] "637.2708_2.5777"                       
    [1874] "639.1931_0.3202"                       
    [1875] "640.1516_2.9411"                       
    [1876] "640.1518_2.8068"                       
    [1877] "641.1508_2.5244"                       
    [1878] "641.151_2.8021"                        
    [1879] "641.1694_2.1918"                       
    [1880] "641.1698_2.3413"                       
    [1881] "642.1546_2.5278"                       
    [1882] "642.1644_0.3096"                       
    [1883] "642.1669_2.9324"                       
    [1884] "643.1431_3.2702"                       
    [1885] "643.1431_3.0357"                       
    [1886] "643.1676_2.1971"                       
    [1887] "645.055_2.7278"                        
    [1888] "645.1424_2.9424"                       
    [1889] "645.1597_3.0613"                       
    [1890] "645.188_1.6162"                        
    [1891] "646.1682_5.8605"                       
    [1892] "646.1914_1.6104"                       
    [1893] "647.0983_0.4005"                       
    [1894] "647.1896_2.4366"                       
    [1895] "647.2039_1.5898"                       
    [1896] "647.2093_3.0611"                       
    [1897] "648.1656_5.8614"                       
    [1898] "648.1839_5.0486"                       
    [1899] "648.184_5.2365"                        
    [1900] "649.1354_2.9307"                       
    [1901] "649.1686_5.8613"                       
    [1902] "649.1825_0.5025"                       
    [1903] "649.2104_2.5797"                       
    [1904] "649.2503_2.8271"                       
    [1905] "649.2672_4.2808"                       
    [1906] "650.1773_2.2374"                       
    [1907] "650.1813_5.2358"                       
    [1908] "651.0545_0.3533"                       
    [1909] "651.1251_2.9745"                       
    [1910] "651.1538_3.0561"                       
    [1911] "651.4108_6.4076"                       
    [1912] "652.9159_0.3661"                       
    [1913] "653.1245_2.7234"                       
    [1914] "653.147_3.0622"                        
    [1915] "655.1119_2.2766"                       
    [1916] "657.0559_2.0373"                       
    [1917] "657.1743_2.5186"                       
    [1918] "658.1563_2.0498"                       
    [1919] "658.1621_2.723"                        
    [1920] "658.1747_2.173"                        
    [1921] "659.2109_2.0499"                       
    [1922] "661.1377_2.3877"                       
    [1923] "661.1538_2.4759"                       
    [1924] "661.1594_2.0493"                       
    [1925] "661.2287_3.6787"                       
    [1926] "663.1014_3.0626"                       
    [1927] "663.1288_2.9313"                       
    [1928] "663.1314_3.176"                        
    [1929] "663.1505_2.7226"                       
    [1930] "663.1983_0.3255"                       
    [1931] "665.1466_2.9748"                       
    [1932] "665.1487_1.6637"                       
    [1933] "665.172_2.6426"                        
    [1934] "665.2127_0.3243"                       
    [1935] "665.3859_2.9516"                       
    [1936] "665.3901_5.8624"                       
    [1937] "667.1404_3.0653"                       
    [1938] "669.0557_0.41"                         
    [1939] "669.0952_2.7226"                       
    [1940] "670.9256_0.4228"                       
    [1941] "671.1035_3.0619"                       
    [1942] "672.1779_3.0616"                       
    [1943] "673.1195_2.0466"                       
    [1944] "673.144_2.6989"                        
    [1945] "674.1513_1.9187"                       
    [1946] "675.131_2.9411"                        
    [1947] "675.2078_2.8606"                       
    [1948] "675.2078_3.1635"                       
    [1949] "676.1631_2.1729"                       
    [1950] "677.1465_2.9329"                       
    [1951] "677.163_2.1729"                        
    [1952] "677.1684_3.0636"                       
    [1953] "677.1777_0.4942"                       
    [1954] "677.2234_3.583"                        
    [1955] "677.2297_2.2213"                       
    [1956] "677.3538_3.9976"                       
    [1957] "678.1711_3.0623"                       
    [1958] "679.1197_2.7223"                       
    [1959] "679.1461_3.0621"                       
    [1960] "679.1684_3.0638"                       
    [1961] "679.1722_1.7959"                       
    [1962] "681.0755_2.7228"                       
    [1963] "681.122_3.1753"                        
    [1964] "681.1642_2.9751"                       
    [1965] "681.1667_2.1218"                       
    [1966] "681.1668_3.3701"                       
    [1967] "681.2394_2.3058"                       
    [1968] "681.2395_2.4598"                       
    [1969] "682.117_2.0508"                        
    [1970] "683.1128_3.0652"                       
    [1971] "683.2247_0.3201"                       
    [1972] "683.2549_2.3258"                       
    [1973] "684.0545_2.7235"                       
    [1974] "685.0611_2.7228"                       
    [1975] "685.1829_1.585"                        
    [1976] "687.0935_0.4011"                       
    [1977] "687.0988_0.6083"                       
    [1978] "688.9353_0.4382"                       
    [1979] "689.0088_0.358"                        
    [1980] "689.1508_2.5286"                       
    [1981] "691.1485_2.9427"                       
    [1982] "693.0634_2.7233"                       
    [1983] "693.0852_2.9966"                       
    [1984] "693.136_3.0656"                        
    [1985] "693.1416_2.7232"                       
    [1986] "693.2765_3.2637"                       
    [1987] "693.3851_3.9283"                       
    [1988] "694.1449_2.7223"                       
    [1989] "695.0609_2.7239"                       
    [1990] "695.091_3.063"                         
    [1991] "695.1377_3.1762"                       
    [1992] "695.1827_2.1997"                       
    [1993] "695.1877_0.3198"                       
    [1994] "695.2078_2.5798"                       
    [1995] "695.2245_0.3199"                       
    [1996] "695.4006_6.0398"                       
    [1997] "695.4007_5.6015"                       
    [1998] "697.129_2.9751"                        
    [1999] "697.2026_0.3116"                       
    [2000] "697.3798_6.5994"                       
    [2001] "697.3799_6.4066"                       
    [2002] "699.0777_3.0653"                       
    [2003] "699.1591_0.4944"                       
    [2004] "703.2662_2.4797"                       
    [2005] "705.042_2.9747"                        
    [2006] "705.1038_0.408"                        
    [2007] "705.1627_0.3256"                       
    [2008] "705.1667_2.0944"                       
    [2009] "705.1668_2.4309"                       
    [2010] "706.4294_6.9725"                       
    [2011] "707.0179_0.4031"                       
    [2012] "707.0788_3.0614"                       
    [2013] "707.101_3.0169"                        
    [2014] "707.1137_2.9423"                       
    [2015] "707.1576_3.0623"                       
    [2016] "707.2398_2.2124"                       
    [2017] "708.0953_3.0634"                       
    [2018] "708.2395_2.2113"                       
    [2019] "708.4612_6.7979"                       
    [2020] "709.0767_3.0611"                       
    [2021] "709.1048_3.0648"                       
    [2022] "709.1292_2.9314"                       
    [2023] "709.1606_2.7213"                       
    [2024] "710.1634_2.7201"                       
    [2025] "711.1105_0.4733"                       
    [2026] "711.2863_2.865"                        
    [2027] "711.2867_3.2593"                       
    [2028] "713.1283_3.1768"                       
    [2029] "715.1333_0.493"                        
    [2030] "715.1335_0.3355"                       
    [2031] "716.1366_0.4927"                       
    [2032] "716.1368_0.3357"                       
    [2033] "717.101_3.0634"                        
    [2034] "717.2361_2.8281"                       
    [2035] "719.1197_0.8202"                       
    [2036] "719.1462_0.5433"                       
    [2037] "719.1503_1.8881"                       
    [2038] "719.2318_2.3256"                       
    [2039] "720.1586_2.2656"                       
    [2040] "720.1587_2.5337"                       
    [2041] "720.4452_7.1809"                       
    [2042] "721.0389_2.7236"                       
    [2043] "722.1702_2.3464"                       
    [2044] "722.1738_2.0921"                       
    [2045] "723.0369_2.723"                        
    [2046] "723.1087_2.3871"                       
    [2047] "723.175_3.0608"                        
    [2048] "723.1928_2.2232"                       
    [2049] "723.2135_3.1509"                       
    [2050] "723.2499_3.1592"                       
    [2051] "724.1961_2.463"                        
    [2052] "725.2658_2.7009"                       
    [2053] "725.2659_2.8832"                       
    [2054] "727.0259_2.7536"                       
    [2055] "727.0881_0.463"                        
    [2056] "727.144_3.1776"                        
    [2057] "727.1935_1.9515"                       
    [2058] "729.0249_2.7542"                       
    [2059] "729.04_2.7219"                         
    [2060] "729.1135_3.0671"                       
    [2061] "729.1644_2.2111"                       
    [2062] "729.2526_4.0738"                       
    [2063] "730.4294_6.6702"                       
    [2064] "731.038_2.7221"                        
    [2065] "731.0965_0.4974"                       
    [2066] "731.1412_2.2108"                       
    [2067] "733.0349_2.7227"                       
    [2068] "733.1407_2.8371"                       
    [2069] "733.1435_0.5068"                       
    [2070] "733.144_0.3299"                        
    [2071] "733.2053_2.8274"                       
    [2072] "735.2002_2.016"                        
    [2073] "735.2033_2.2556"                       
    [2074] "735.2035_2.0806"                       
    [2075] "735.2199_1.5913"                       
    [2076] "735.2864_4.193"                        
    [2077] "736.1596_2.6606"                       
    [2078] "736.2232_1.5905"                       
    [2079] "737.1242_3.0332"                       
    [2080] "737.1714_1.9431"                       
    [2081] "737.1718_1.1034"                       
    [2082] "737.1733_2.9536"                       
    [2083] "737.1928_2.1905"                       
    [2084] "737.1928_2.8923"                       
    [2085] "737.193_2.5524"                        
    [2086] "737.2354_1.5594"                       
    [2087] "738.1959_2.1896"                       
    [2088] "739.1336_2.5391"                       
    [2089] "739.1399_3.0648"                       
    [2090] "739.1665_2.6864"                       
    [2091] "739.1872_1.7731"                       
    [2092] "739.2055_2.5241"                       
    [2093] "741.1672_2.7203"                       
    [2094] "741.2969_2.9927"                       
    [2095] "742.0674_0.6041"                       
    [2096] "743.0488_3.0623"                       
    [2097] "745.0462_3.062"                        
    [2098] "745.0988_0.417"                        
    [2099] "745.1357_2.2108"                       
    [2100] "746.1022_0.4164"                       
    [2101] "747.0449_3.0619"                       
    [2102] "747.1403_2.2109"                       
    [2103] "747.2864_3.1853"                       
    [2104] "748.0999_0.5836"                       
    [2105] "749.1063_0.5124"                       
    [2106] "749.2115_0.3186"                       
    [2107] "753.1192_2.4923"                       
    [2108] "753.2243_3.098"                        
    [2109] "755.1812_3.0612"                       
    [2110] "755.2034_2.525"                        
    [2111] "757.2191_2.4345"                       
    [2112] "759.1692_2.2279"                       
    [2113] "761.2623_3.2637"                       
    [2114] "763.1093_0.5229"                       
    [2115] "763.1899_0.312"                        
    [2116] "763.2147_1.5619"                       
    [2117] "765.0498_0.355"                        
    [2118] "765.1778_1.8244"                       
    [2119] "765.1779_2.2381"                       
    [2120] "765.2297_1.5512"                       
    [2121] "766.1613_1.879"                        
    [2122] "766.1621_2.296"                        
    [2123] "766.2256_0.3255"                       
    [2124] "768.4662_6.2658"                       
    [2125] "769.1826_2.306"                        
    [2126] "769.1829_2.6175"                       
    [2127] "769.2191_2.6218"                       
    [2128] "771.1013_2.2112"                       
    [2129] "771.198_2.3088"                        
    [2130] "771.1985_2.0653"                       
    [2131] "772.1873_1.5923"                       
    [2132] "773.172_2.7224"                        
    [2133] "776.1985_2.1657"                       
    [2134] "776.7_2.1657"                          
    [2135] "777.09_3.0696"                         
    [2136] "777.2305_1.6384"                       
    [2137] "777.2317_3.2636"                       
    [2138] "778.2338_1.639"                        
    [2139] "779.1939_2.2608"                       
    [2140] "779.2077_2.8271"                       
    [2141] "779.2457_1.6426"                       
    [2142] "779.3028_1.6441"                       
    [2143] "780.2492_1.6428"                       
    [2144] "783.0563_0.4017"                       
    [2145] "783.0683_2.044"                        
    [2146] "785.2127_2.3663"                       
    [2147] "785.2136_2.1376"                       
    [2148] "787.1875_3.0627"                       
    [2149] "787.193_2.3323"                        
    [2150] "787.1936_1.9952"                       
    [2151] "791.1626_2.7227"                       
    [2152] "791.2667_2.943"                        
    [2153] "793.2405_2.2059"                       
    [2154] "793.4013_5.3499"                       
    [2155] "795.3077_3.9454"                       
    [2156] "795.4893_6.4576"                       
    [2157] "795.4895_6.1824"                       
    [2158] "797.0785_2.2104"                       
    [2159] "799.4244_6.7949"                       
    [2160] "801.0647_0.42"                         
    [2161] "803.0088_0.3613"                       
    [2162] "803.1856_2.7235"                       
    [2163] "805.1779_3.0634"                       
    [2164] "805.2041_2.1777"                       
    [2165] "807.1749_2.0634"                       
    [2166] "809.1531_2.7229"                       
    [2167] "811.3025_2.97"                         
    [2168] "813.105_0.4928"                        
    [2169] "813.2358_0.3245"                       
    [2170] "813.2453_2.1566"                       
    [2171] "817.183_2.1778"                        
    [2172] "817.1998_2.0417"                       
    [2173] "820.4721_4.6174"                       
    [2174] "821.152_0.3525"                        
    [2175] "823.1687_2.7229"                       
    [2176] "823.1702_1.9971"                       
    [2177] "823.2338_3.2636"                       
    [2178] "825.1698_1.9968"                       
    [2179] "827.18_2.9751"                         
    [2180] "834.215_2.285"                         
    [2181] "837.1644_2.9419"                       
    [2182] "837.1838_3.0629"                       
    [2183] "837.1842_2.7232"                       
    [2184] "839.2304_0.4709"                       
    [2185] "839.4065_3.3307"                       
    [2186] "851.1997_3.0615"                       
    [2187] "855.1748_3.0635"                       
    [2188] "855.175_2.7249"                        
    [2189] "861.1668_2.2754"                       
    [2190] "863.1997_1.643"                        
    [2191] "864.6919_2.4025"                       
    [2192] "864.692_2.5795"                        
    [2193] "Procyanidin C1"                        
    [2194] "866.2009_2.7398"                       
    [2195] "866.2015_2.2923"                       
    [2196] "867.2619_1.5665"                       
    [2197] "869.1905_3.0624"                       
    [2198] "869.2775_1.5752"                       
    [2199] "879.1769_2.0725"                       
    [2200] "879.1775_2.3634"                       
    [2201] "879.2126_2.6994"                       
    [2202] "879.2132_2.9716"                       
    [2203] "882.2019_2.0642"                       
    [2204] "883.1652_2.5891"                       
    [2205] "885.2663_2.5574"                       
    [2206] "889.1427_0.4909"                       
    [2207] "895.1918_0.3516"                       
    [2208] "895.1935_2.1839"                       
    [2209] "899.2438_2.4536"                       
    [2210] "899.245_2.1284"                        
    [2211] "900.2359_2.4142"                       
    [2212] "901.2385_2.4168"                       
    [2213] "903.1518_2.4615"                       
    [2214] "907.152_0.5042"                        
    [2215] "909.1847_0.5792"                       
    [2216] "911.1283_0.4906"                       
    [2217] "911.2876_1.659"                        
    [2218] "925.1626_0.5839"                       
    [2219] "926.1659_0.5836"                       
    [2220] "927.102_0.4895"                        
    [2221] "927.1676_0.5841"                       
    [2222] "929.2139_2.5395"                       
    [2223] "929.2709_3.4011"                       
    [2224] "933.2302_2.2039"                       
    [2225] "933.2494_2.2969"                       
    [2226] "933.2512_2.0707"                       
    [2227] "934.0712_2.5127"                       
    [2228] "941.1349_0.584"                        
    [2229] "943.0686_0.4893"                       
    [2230] "946.2187_2.0527"                       
    [2231] "947.1474_0.5668"                       
    [2232] "947.2457_2.6649"                       
    [2233] "947.2459_2.348"                        
    [2234] "949.2252_2.1205"                       
    [2235] "949.2463_2.0095"                       
    [2236] "951.5675_7.45"                         
    [2237] "955.0782_0.3511"                       
    [2238] "963.1209_0.5603"                       
    [2239] "963.2407_2.2613"                       
    [2240] "963.2767_2.4711"                       
    [2241] "964.1237_0.5588"                       
    [2242] "965.2499_2.2547"                       
    [2243] "965.2564_2.5613"                       
    [2244] "966.2598_2.5609"                       
    [2245] "967.5625_7.0719"                       
    [2246] "969.2935_1.6011"                       
    [2247] "977.2565_2.3172"                       
    [2248] "978.0796_0.5839"                       
    [2249] "979.0871_0.5826"                       
    [2250] "979.2362_2.1013"                       
    [2251] "981.2032_0.336"                        
    [2252] "981.2511_2.8825"                       
    [2253] "993.2515_2.2426"                       
    [2254] "1007.331_0.3252"                       
    [2255] "1008.7234_2.6265"                      
    [2256] "1009.7263_2.3996"                      
    [2257] "1010.228_2.6269"                       
    [2258] "1027.2295_2.8053"                      
    [2259] "1028.3691_3.0675"                      
    [2260] "1029.228_2.2432"                       
    [2261] "1095.2829_2.1867"                      
    [2262] "1109.2985_2.3128"                      
    [2263] "1111.2778_2.1126"                      
    [2264] "1125.2934_2.2365"                      
    [2265] "1139.3091_2.2804"                      
    [2266] "1143.2574_0.3316"                      
    [2267] "1151.2453_1.9834"                      
    [2268] "1152.7546_2.4391"                      
    [2269] "1153.2612_2.5123"                      
    [2270] "1153.302_2.6239"                       
    [2271] "1153.3036_2.8361"                      
    [2272] "1155.3033_2.2112"                      
    [2273] "1171.3139_2.8644"                      
    [2274] "1183.3138_2.8513"                      
    [2275] "1201.3241_2.8854"                      
    [2276] "1241.3522_3.0614"                      
    [2277] "1243.3452_3.0614"                      
    [2278] "1257.3212_3.0619"                      
    [2279] "1272.2832_3.062"                       
    [2280] "1273.2868_3.062"                       
    [2281] "1287.3457_2.2017"                      
    [2282] "1293.2954_3.0623"                      
    [2283] "1296.7852_2.4929"                      
    [2284] "1333.3661_2.5313"                      
    [2285] "1509.3982_2.2087"                      

``` r
# find which features have no variance across fresh juice and MEF
# having samples with no variance won't work for a heatmap
all_MEF_comparisons_novariance <- feature_table_long_log2_anova %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# what features have no variance?
all_MEF_comparisons_novariance$id
```

    character(0)

``` r
# wrangle

feature_table_wide_log2_anova <- feature_table_long_log2_anova %>%
  select(-"id")  %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")


feature_table_MEF_anova_for_heatmap <- feature_table_wide_log2_anova %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  select(-sample_or_qc) %>% # remove columns we don't need
  column_to_rownames(var = "sample_name") # move the sample names to rownames for plotting


heatmap <- pheatmap(feature_table_MEF_anova_for_heatmap[,-1], # remove treatment column
         cluster_rows = TRUE, # HCA on samples
         cluster_cols = TRUE, # HCA on features 
         scale = "column",
         show_colnames = FALSE, # scale for each feature
         legend_breaks = c(-3, 0, 3),
         legend_labels = c('Low relative abundance', 'Medium relative abundance', 'High relative abundance'),
         border_color = F,
         filename = "heatmap_anova.png")
```

``` r
ncol(feature_table_MEF_for_heatmap)
```

    [1] 2283

``` r
ncol(feature_table_MEF_anova_for_heatmap)
```

    [1] 1714
