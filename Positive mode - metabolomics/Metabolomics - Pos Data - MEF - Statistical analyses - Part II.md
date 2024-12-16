# Metabolomics - Pos Data - MEF - Statistical analyses - Part II
Giovana Domeneghini Mercali

# Load libraries

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
library(rcompanion)
library(purrr)
library(multcompView) # Load the multcompView package for comparison letters

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

# Read in data

Positive mode.

``` r
# read in metabolomics data

feature_table_long_log2 <- read_csv("feature_table_long_log2__annotated_pos.csv",
                  trim_ws = TRUE) 
```

    Rows: 81411 Columns: 10
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (6): sample_name, treatment, sample_or_qc, mz_rt, Compound, id
    dbl (4): run_order, mz, rt, rel_abund

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data
metadata <- read_csv("metadata_pos_mef.csv",
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
annotation <- read_csv("Annotation_pos.csv", 
                  locale=locale(encoding = "latin1"))
```

    Rows: 24 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (2): Compound, mz_rt

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Compariosn between compounds D

## Boxplots

Filter to select the compound you want to plot.

``` r
compound_D <- feature_table_long_log2 %>%
  filter(id == "Cyanidin hexoside" | id == "Cyanidin malonyl hexoside" | id == "Cyanidin rutinoside" | id == "N-Caffeoylputrescine" | id == "p-Coumaroylputrescine" | id == "Phenylalanine") %>%
   filter(sample_or_qc == "Sample")
```

Then plot.

``` r
boxplot_compound_D <- compound_D %>%
  ggplot(aes(x = factor(treatment, levels = c("FJ", "MT", "MEF", "MEF+SS")), 
             y = rel_abund, fill = treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_wrap(vars(id), nrow = 2, scales = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5")) +
  theme(legend.position = "none") +
  labs(x = "Treatment",
  y = "log2(relative abundance)",
  title = "")

boxplot_compound_D
```

![](Metabolomics---Pos-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-4-1.png)

## Testing assumptions

### Normality

``` r
# testing normality by group
compound_D %>%
  drop_na(rel_abund) %>% # remove NAs
  group_by(id) %>% # test by species
  shapiro_test(rel_abund) # test for normality
```

    # A tibble: 6 × 4
      id                        variable  statistic        p
      <chr>                     <chr>         <dbl>    <dbl>
    1 Cyanidin hexoside         rel_abund     0.851 0.00442 
    2 Cyanidin malonyl hexoside rel_abund     0.890 0.0222  
    3 Cyanidin rutinoside       rel_abund     0.837 0.00260 
    4 N-Caffeoylputrescine      rel_abund     0.759 0.000165
    5 Phenylalanine             rel_abund     0.766 0.000208
    6 p-Coumaroylputrescine     rel_abund     0.765 0.000202

### Constant variance

``` r
# testing constant variance

compound_D %>%
  drop_na(rel_abund) %>% # remove NAs
  levene_test(rel_abund ~ id) # test for constant variance
```

    Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    factor.

    # A tibble: 1 × 4
        df1   df2 statistic        p
      <int> <int>     <dbl>    <dbl>
    1     5   120      4.97 0.000360

``` r
compound_D %>%
  drop_na(rel_abund) %>%
  ggplot(aes(x = rel_abund, y = id, fill = id)) +
  ggridges::geom_density_ridges(alpha = 0.7) +
  # scale_y_discrete(labels = c("Cyanidin hexoside", "Cyanidin malonyl hexoside", "Cyanidin rutinoside", "N-Caffeoylputrescine", "p-Coumaroylputrescine", "Phenylalanine")) +
  theme_minimal() +  
  theme(legend.position = "none") +
  labs(x = "Relative abundance",
       y = "Compound",
       title = "Distribution of relative abundace for some compounds in blended juice")
```

    Picking joint bandwidth of 0.702

![](Metabolomics---Pos-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-7-1.png)

Find which features have no variance across treatments

``` r
treatments_novariance <- compound_D %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)
```

No features with zero variance.

## Nonparametric tests

### kruskal Wallis test

``` r
kruskal_compound_D <- compound_D %>%
  drop_na(rel_abund, treatment) %>%            # Remove NA values first
  group_by(id) %>%                             # Group by 'id'
  kruskal_test(rel_abund ~ treatment) %>%      # Apply the Kruskal-Wallis test
  add_significance() %>%     
  arrange(p)     
           
kruskal_compound_D
```

    # A tibble: 6 × 8
      id                        .y.        n statistic    df       p method p.signif
      <chr>                     <chr>  <int>     <dbl> <int>   <dbl> <chr>  <chr>   
    1 Cyanidin hexoside         rel_a…    21      18.8     3 3.06e-4 Krusk… ***     
    2 Cyanidin malonyl hexoside rel_a…    21      18.8     3 3.06e-4 Krusk… ***     
    3 Cyanidin rutinoside       rel_a…    21      18.8     3 3.06e-4 Krusk… ***     
    4 p-Coumaroylputrescine     rel_a…    21      17.7     3 5.01e-4 Krusk… ***     
    5 N-Caffeoylputrescine      rel_a…    21      17.2     3 6.45e-4 Krusk… ***     
    6 Phenylalanine             rel_a…    21      15.2     3 1.62e-3 Krusk… **      

### Post-hoc test

Run Dunn test

``` r
(kruskal_compound_D_posthoc <- compound_D %>%
      drop_na(rel_abund, treatment) %>%        # Remove NA values first
      group_by(id) %>%                         # Group by 'id'
      dunn_test(rel_abund ~ treatment,
                p.adjust.method = "BH"))
```

    # A tibble: 36 × 10
       id     .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
     * <chr>  <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
     1 Cyani… rel_… FJ     MEF        6     5     -2.79 5.20e-3 1.56e-2 *           
     2 Cyani… rel_… FJ     MEF+SS     6     5     -4.13 3.70e-5 2.22e-4 ***         
     3 Cyani… rel_… FJ     MT         6     5     -1.46 1.43e-1 2.03e-1 ns          
     4 Cyani… rel_… MEF    MEF+SS     5     5     -1.27 2.03e-1 2.03e-1 ns          
     5 Cyani… rel_… MEF    MT         5     5      1.27 2.03e-1 2.03e-1 ns          
     6 Cyani… rel_… MEF+SS MT         5     5      2.55 1.08e-2 2.17e-2 *           
     7 Cyani… rel_… FJ     MEF        6     5     -2.79 5.20e-3 1.56e-2 *           
     8 Cyani… rel_… FJ     MEF+SS     6     5     -4.13 3.70e-5 2.22e-4 ***         
     9 Cyani… rel_… FJ     MT         6     5     -1.46 1.43e-1 2.03e-1 ns          
    10 Cyani… rel_… MEF    MEF+SS     5     5     -1.27 2.03e-1 2.03e-1 ns          
    # ℹ 26 more rows

``` r
# Create a column called comparison

kruskal_compound_D_posthocc_1 <- kruskal_compound_D_posthoc %>%
  mutate(comparison = glue("{group1} - {group2}")) 
  
# Group data by 'id' and create a nested data frame
group_cldList <- kruskal_compound_D_posthocc_1 %>%
  group_by(id) %>%
  nest() %>%  # Nest data for each 'id'
  
# Apply cldList function on each nested data
mutate(cld_results = map(data, ~ cldList(p.adj ~ comparison, data = .x, threshold = 0.05))) %>%
 unnest(cld_results) %>%  # Unnest the results
  select(-data) %>%       # Remove the nested column and ungroup the result
  ungroup()

# Display the result
group_cldList <- group_cldList %>%
  rename(treatment = Group) 
```

Make a table that has the maximum relative abundance for each treatment.

``` r
# Calculate maximum rel_abund for each id and treatment
max_relabund <- compound_D %>%
  drop_na(rel_abund, treatment) %>%  # Remove NA values
  group_by(id, treatment) %>%        # Group by 'id' and 'treatment'
  summarize(max_tot_relabund = max(rel_abund, na.rm = TRUE)) %>%  # Calculate max
  ungroup()  # Remove grouping
```

    `summarise()` has grouped output by 'id'. You can override using the `.groups`
    argument.

``` r
# Display the result
max_relabund
```

    # A tibble: 24 × 3
       id                        treatment max_tot_relabund
       <chr>                     <chr>                <dbl>
     1 Cyanidin hexoside         FJ                    24.5
     2 Cyanidin hexoside         MEF                   23.5
     3 Cyanidin hexoside         MEF+SS                22.3
     4 Cyanidin hexoside         MT                    24.1
     5 Cyanidin malonyl hexoside FJ                    20.1
     6 Cyanidin malonyl hexoside MEF                   18.3
     7 Cyanidin malonyl hexoside MEF+SS                17.0
     8 Cyanidin malonyl hexoside MT                    19.1
     9 Cyanidin rutinoside       FJ                    21.6
    10 Cyanidin rutinoside       MEF                   21.1
    # ℹ 14 more rows

Bind the groups to the maximum relative abundance.

``` r
dunn_compound_D_for_plotting <- max_relabund %>%
  left_join(group_cldList, by = c("id", "treatment"))
```

Plot

``` r
boxplot_dunn <- boxplot_compound_D +
  geom_text(data = dunn_compound_D_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = dunn_compound_D_for_plotting$Letter)) +
  labs(caption = "For each compound, treatmnets with different letters are significant different unsing Krushal Wallis test")

boxplot_dunn 
```

    Warning: Use of `dunn_compound_D_for_plotting$Letter` is discouraged.
    ℹ Use `Letter` instead.

![](Metabolomics---Pos-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-14-1.png)

## Parametric tests

### ANOVA test

``` r
# Prepare your data 
clean_data <- compound_D %>%
  drop_na(rel_abund, treatment) %>%
  group_by(id) %>%
  nest()  # Nest data for each 'id'

# Apply ANOVA to each nested data frame
results <- clean_data %>%
  mutate(
    anova_model = map(data, ~ aov(rel_abund ~ treatment, data = .x)),
    anova_summary = map(anova_model, tidy)  # Use broom::tidy to summarize ANOVA results
  )

# Print results for each id
results %>%
  select(id, anova_summary) %>%
  unnest(anova_summary) %>%
  print()
```

    # A tibble: 12 × 7
    # Groups:   id [6]
       id                        term         df   sumsq  meansq statistic   p.value
       <chr>                     <chr>     <dbl>   <dbl>   <dbl>     <dbl>     <dbl>
     1 Phenylalanine             treatment     3  29.7    9.91       117.   1.51e-11
     2 Phenylalanine             Residuals    17   1.44   0.0848      NA   NA       
     3 p-Coumaroylputrescine     treatment     3  91.2   30.4        126.   8.24e-12
     4 p-Coumaroylputrescine     Residuals    17   4.10   0.241       NA   NA       
     5 N-Caffeoylputrescine      treatment     3 106.    35.5        135.   4.76e-12
     6 N-Caffeoylputrescine      Residuals    17   4.48   0.263       NA   NA       
     7 Cyanidin hexoside         treatment     3  16.4    5.45        98.0  6.25e-11
     8 Cyanidin hexoside         Residuals    17   0.946  0.0556      NA   NA       
     9 Cyanidin malonyl hexoside treatment     3  32.4   10.8        180.   4.36e-13
    10 Cyanidin malonyl hexoside Residuals    17   1.02   0.0598      NA   NA       
    11 Cyanidin rutinoside       treatment     3   6.73   2.24       118.   1.38e-11
    12 Cyanidin rutinoside       Residuals    17   0.323  0.0190      NA   NA       

``` r
# To perform Tukey HSD test and print results
results <- results %>%
  mutate(
    tukey_test = map(anova_model, TukeyHSD),
    tukey_summary = map(tukey_test, tidy)  # Use broom::tidy to summarize Tukey HSD results
  )

# Print Tukey HSD results
table_results <- results %>%
  select(id, tukey_summary) %>%
  unnest(tukey_summary) %>%
  arrange(adj.p.value)%>%
  print()
```

    # A tibble: 36 × 8
    # Groups:   id [6]
       id          term  contrast null.value estimate conf.low conf.high adj.p.value
       <chr>       <chr> <chr>         <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
     1 Cyanidin m… trea… MEF+SS-…          0    -3.28   -3.70      -2.86    5.53e-13
     2 Cyanidin r… trea… MEF+SS-…          0    -1.50   -1.74      -1.27    8.67e-12
     3 Cyanidin h… trea… MEF+SS-…          0    -2.29   -2.69      -1.88    6.18e-11
     4 Cyanidin m… trea… MT-MEF+…          0     2.36    1.92       2.80    1.37e-10
     5 N-Caffeoyl… trea… MEF+SS-…          0    -4.71   -5.59      -3.82    1.52e-10
     6 p-Coumaroy… trea… MEF+SS-…          0    -4.42   -5.27      -3.58    2.05e-10
     7 N-Caffeoyl… trea… MT-MEF+…          0     4.62    3.70       5.54    4.09e-10
     8 N-Caffeoyl… trea… MEF-FJ            0    -4.38   -5.26      -3.50    4.77e-10
     9 Cyanidin r… trea… MT-MEF+…          0     1.22    0.972      1.47    5.30e-10
    10 p-Coumaroy… trea… MT-MEF+…          0     4.30    3.42       5.19    6.21e-10
    # ℹ 26 more rows

``` r
# Create a column called comparison

table_results_1 <- table_results %>%
rename(comparison = `contrast`)

# Group data by 'id' and create a nested data frame
group_cldList_2 <- table_results_1 %>%
  group_by(id) %>%
  nest() %>%  # Nest data for each 'id'
  
# Apply cldList function on each nested data
mutate(cld_results = map(data, ~ cldList(adj.p.value ~ comparison, data = .x, threshold = 0.05))) %>%
 unnest(cld_results) %>%  # Unnest the results
  select(-data) %>%       # Remove the nested column and ungroup the result
  ungroup()


# Display the result
group_cldList_2 <- group_cldList_2 %>%
  rename(treatment = Group) 
```

Join the groups to the maximum relative abundance.

``` r
tukey_compound_D_for_plotting <- max_relabund %>%
  left_join(group_cldList_2, by = c("id", "treatment"))
```

``` r
Boxplot_tukey_compound_D <- boxplot_compound_D +
  geom_text(data = tukey_compound_D_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = Letter)) +
  labs(title = "", tag = "C")

Boxplot_tukey_compound_D
```

![](Metabolomics---Pos-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
ggsave(plot = Boxplot_tukey_compound_D,
       filename = "boxplot_compound_D_pos.svg")
```

    Saving 7 x 5 in image
