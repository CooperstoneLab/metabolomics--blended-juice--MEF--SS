# Metabolomics - Neg Data - MEF - Statistical analyses - Part II
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

Negative mode.

``` r
# read in metabolomics data

feature_table_long_log2 <- read_csv("feature_table_long_log2__annotated_neg.csv",
                  trim_ws = TRUE) 
```

    Rows: 75405 Columns: 10
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (6): sample_name, treatment, sample_or_qc, mz_rt, Compound, id
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

# Compariosn between compounds A

## Boxplots

Filter to select the compound you want to plot.

``` r
compound_A <- feature_table_long_log2 %>%
  filter(id == "3,4-Dihydroxybenzoic acid" | id == "Gallic acid" | id == "Ferulic acid" | id == "Caffeic acid" | id =="2,4,6-Trihydroxybenzoic acid" | id == "Ellagic acid") %>%
   filter(sample_or_qc == "Sample")
```

Then plot.

``` r
boxplot_compound_A <- compound_A %>%
  ggplot(aes(x = factor(treatment, levels = c("FJ", "MT", "MEF", "MEF+SS")), 
             y = rel_abund, fill = treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_wrap(vars(id), nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5")) +
  labs(x = "Treatment",
  y = "log2(relative abundance)",
  title = "")

boxplot_compound_A
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-4-1.png)

## Testing assumptions

### Normality

``` r
# testing normality by group
compound_A %>%
  drop_na(rel_abund) %>% # remove NAs
  group_by(id) %>% # test by species
  shapiro_test(rel_abund) # test for normality
```

    # A tibble: 6 × 4
      id                           variable  statistic       p
      <chr>                        <chr>         <dbl>   <dbl>
    1 2,4,6-Trihydroxybenzoic acid rel_abund     0.848 0.00386
    2 3,4-Dihydroxybenzoic acid    rel_abund     0.878 0.0136 
    3 Caffeic acid                 rel_abund     0.886 0.0190 
    4 Ellagic acid                 rel_abund     0.859 0.00621
    5 Ferulic acid                 rel_abund     0.880 0.0148 
    6 Gallic acid                  rel_abund     0.859 0.00620

### Constant variance

``` r
# testing constant variance

compound_A %>%
  drop_na(rel_abund) %>% # remove NAs
  levene_test(rel_abund ~ id) # test for constant variance
```

    Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    factor.

    # A tibble: 1 × 4
        df1   df2 statistic          p
      <int> <int>     <dbl>      <dbl>
    1     5   120      7.38 0.00000455

``` r
compound_A %>%
  drop_na(rel_abund) %>%
  ggplot(aes(x = rel_abund, y = id, fill = id)) +
  ggridges::geom_density_ridges(alpha = 0.7) +
  # scale_y_discrete(labels = c("3,4-Dihydroxybenzoic acid", "Gallic acid", "Ferulic acid", "Caffeic acid", "2,4,6-Trihydroxybenzoic acid", "Ellagic acid")) +
  theme_minimal() +  
  theme(legend.position = "none") +
  labs(x = "Relative abundance",
       y = "Compound",
       title = "Distribution of relative abundace for some compounds in blended juice")
```

    Picking joint bandwidth of 0.611

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-7-1.png)

Find which features have no variance across treatments

``` r
treatments_novariance <- compound_A %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)
```

No features with zero variance.

## Nonparametric tests

### kruskal Wallis test

``` r
kruskal_compound_A <- compound_A %>%
  drop_na(rel_abund, treatment) %>%            # Remove NA values first
  group_by(id) %>%                             # Group by 'id'
  kruskal_test(rel_abund ~ treatment) %>%      # Apply the Kruskal-Wallis test
  add_significance() %>%     
  arrange(p)     
           
kruskal_compound_A
```

    # A tibble: 6 × 8
      id                         .y.       n statistic    df       p method p.signif
      <chr>                      <chr> <int>     <dbl> <int>   <dbl> <chr>  <chr>   
    1 3,4-Dihydroxybenzoic acid  rel_…    21      18.5     3 3.44e-4 Krusk… ***     
    2 Caffeic acid               rel_…    21      17.5     3 5.54e-4 Krusk… ***     
    3 Ellagic acid               rel_…    21      17.2     3 6.53e-4 Krusk… ***     
    4 2,4,6-Trihydroxybenzoic a… rel_…    21      17.1     3 6.77e-4 Krusk… ***     
    5 Ferulic acid               rel_…    21      16.8     3 7.67e-4 Krusk… ***     
    6 Gallic acid                rel_…    21      16.4     3 9.46e-4 Krusk… ***     

### Post-hoc test

Run Dunn test

``` r
(kruskal_compound_A_posthoc <- compound_A %>%
      drop_na(rel_abund, treatment) %>%        # Remove NA values first
      group_by(id) %>%                         # Group by 'id'
      dunn_test(rel_abund ~ treatment,
                p.adjust.method = "BH"))
```

    # A tibble: 36 × 10
       id     .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
     * <chr>  <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
     1 2,4,6… rel_… FJ     MEF        6     5     3.15  1.64e-3 4.91e-3 **          
     2 2,4,6… rel_… FJ     MEF+SS     6     5     1.82  6.90e-2 1.03e-1 ns          
     3 2,4,6… rel_… FJ     MT         6     5    -0.683 4.95e-1 4.95e-1 ns          
     4 2,4,6… rel_… MEF    MEF+SS     5     5    -1.27  2.03e-1 2.43e-1 ns          
     5 2,4,6… rel_… MEF    MT         5     5    -3.67  2.43e-4 1.46e-3 **          
     6 2,4,6… rel_… MEF+SS MT         5     5    -2.40  1.66e-2 3.32e-2 *           
     7 3,4-D… rel_… FJ     MEF        6     5     4.07  4.66e-5 2.79e-4 ***         
     8 3,4-D… rel_… FJ     MEF+SS     6     5     2.85  4.40e-3 1.32e-2 *           
     9 3,4-D… rel_… FJ     MT         6     5     1.46  1.43e-1 2.15e-1 ns          
    10 3,4-D… rel_… MEF    MEF+SS     5     5    -1.17  2.41e-1 2.41e-1 ns          
    # ℹ 26 more rows

``` r
# Create a column called comparison

kruskal_compound_A_posthocc_1 <- kruskal_compound_A_posthoc %>%
  mutate(comparison = glue("{group1} - {group2}")) 
  
# Group data by 'id' and create a nested data frame
group_cldList <- kruskal_compound_A_posthocc_1 %>%
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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Make a table that has the maximum relative abundance for each treatment.

``` r
# Calculate maximum rel_abund for each id and treatment
max_relabund <- compound_A %>%
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
       id                           treatment max_tot_relabund
       <chr>                        <chr>                <dbl>
     1 2,4,6-Trihydroxybenzoic acid FJ                    14.2
     2 2,4,6-Trihydroxybenzoic acid MEF                   18.7
     3 2,4,6-Trihydroxybenzoic acid MEF+SS                17.7
     4 2,4,6-Trihydroxybenzoic acid MT                    14.0
     5 3,4-Dihydroxybenzoic acid    FJ                    18.1
     6 3,4-Dihydroxybenzoic acid    MEF                   21.4
     7 3,4-Dihydroxybenzoic acid    MEF+SS                20.6
     8 3,4-Dihydroxybenzoic acid    MT                    19.1
     9 Caffeic acid                 FJ                    15.9
    10 Caffeic acid                 MEF                   16.8
    # ℹ 14 more rows

Bind the groups to the maximum relative abundance.

``` r
dunn_compound_A_for_plotting <- max_relabund %>%
  left_join(group_cldList, by = c("id", "treatment"))
```

Plot

``` r
boxplot_dunn <- boxplot_compound_A +
  geom_text(data = dunn_compound_A_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = dunn_compound_A_for_plotting$Letter)) +
  labs(caption = "For each compound, treatmnets with different letters are significant different unsing Krushal Wallis test")

boxplot_dunn 
```

    Warning: Use of `dunn_compound_A_for_plotting$Letter` is discouraged.
    ℹ Use `Letter` instead.

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-14-1.png)

## Parametric tests

### ANOVA test

``` r
# Prepare your data 
clean_data <- compound_A %>%
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
       id                           term       df  sumsq  meansq statistic   p.value
       <chr>                        <chr>   <dbl>  <dbl>   <dbl>     <dbl>     <dbl>
     1 3,4-Dihydroxybenzoic acid    treatm…     3 29.9    9.96       304.   5.86e-15
     2 3,4-Dihydroxybenzoic acid    Residu…    17  0.558  0.0328      NA   NA       
     3 2,4,6-Trihydroxybenzoic acid treatm…     3 99.4   33.1        163.   9.83e-13
     4 2,4,6-Trihydroxybenzoic acid Residu…    17  3.45   0.203       NA   NA       
     5 Gallic acid                  treatm…     3 38.2   12.7         95.7  7.57e-11
     6 Gallic acid                  Residu…    17  2.26   0.133       NA   NA       
     7 Caffeic acid                 treatm…     3  7.51   2.50        50.5  1.12e- 8
     8 Caffeic acid                 Residu…    17  0.842  0.0496      NA   NA       
     9 Ferulic acid                 treatm…     3 10.8    3.59        72.7  6.64e-10
    10 Ferulic acid                 Residu…    17  0.839  0.0494      NA   NA       
    11 Ellagic acid                 treatm…     3 25.5    8.50       126.   8.17e-12
    12 Ellagic acid                 Residu…    17  1.15   0.0674      NA   NA       

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
     1 3,4-Dihydr… trea… MEF-FJ            0     2.94     2.63      3.26    1.76e-13
     2 3,4-Dihydr… trea… MEF+SS-…          0     2.49     2.18      2.80    4.27e-13
     3 2,4,6-Trih… trea… MT-MEF            0    -5.15    -5.96     -4.34    8.24e-12
     4 2,4,6-Trih… trea… MEF-FJ            0     4.75     3.97      5.52    1.52e-11
     5 3,4-Dihydr… trea… MT-MEF            0    -1.90    -2.23     -1.58    3.31e-11
     6 Ellagic ac… trea… MT-MEF+…          0     2.57     2.10      3.03    9.00e-11
     7 Ellagic ac… trea… MT-MEF            0     2.47     2.00      2.93    1.74e-10
     8 Gallic acid trea… MEF-FJ            0     3.05     2.43      3.68    6.41e-10
     9 2,4,6-Trih… trea… MT-MEF+…          0    -3.77    -4.58     -2.96    1.25e- 9
    10 Ferulic ac… trea… MT-MEF            0     1.84     1.44      2.24    1.50e- 9
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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Join the groups to the maximum relative abundance.

``` r
tukey_compound_A_for_plotting <- max_relabund %>%
  left_join(group_cldList_2, by = c("id", "treatment"))
```

``` r
Boxplot_tukey_compound_A <- boxplot_compound_A +
  geom_text(data = tukey_compound_A_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = Letter)) +
  labs(title = "", tag = "D")

Boxplot_tukey_compound_A
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
ggsave(plot = Boxplot_tukey_compound_A,
       filename = "boxplot_compound_A_neg.svg")
```

    Saving 7 x 5 in image

# Compariosn between compounds B

## Boxplots

Filter to select the compound you want to plot.

``` r
compound_B <- feature_table_long_log2 %>%
  filter(id == "Epicatechin" |  id == "Procyanidin B2" | id == "Procyanidin C1" | id == "Phloretin" |id =="Phlorhizin" | id == "Quercitrin") %>%
   filter(sample_or_qc == "Sample")
```

Then plot.

``` r
boxplot_compound_B <- compound_B %>%
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
```

## Testing assumptions

### Normality

``` r
# testing normality by group
compound_B %>%
  drop_na(rel_abund) %>% # remove NAs
  group_by(id) %>% # test by species
  shapiro_test(rel_abund) # test for normality
```

    # A tibble: 6 × 4
      id             variable  statistic        p
      <chr>          <chr>         <dbl>    <dbl>
    1 Epicatechin    rel_abund     0.811 0.000985
    2 Phloretin      rel_abund     0.920 0.0865  
    3 Phlorhizin     rel_abund     0.830 0.00198 
    4 Procyanidin B2 rel_abund     0.877 0.0131  
    5 Procyanidin C1 rel_abund     0.887 0.0196  
    6 Quercitrin     rel_abund     0.838 0.00261 

### Constant variance

``` r
# testing constant variance

compound_B %>%
  drop_na(rel_abund) %>% # remove NAs
  levene_test(rel_abund ~ id) # test for constant variance
```

    Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    factor.

    # A tibble: 1 × 4
        df1   df2 statistic         p
      <int> <int>     <dbl>     <dbl>
    1     5   120      6.63 0.0000175

``` r
compound_B %>%
  drop_na(rel_abund) %>%
  ggplot(aes(x = rel_abund, y = id, fill = id)) +
  ggridges::geom_density_ridges(alpha = 0.7) +
  # scale_y_discrete(labels = c("Epicatechin", "Procyanidin B2", "Procyanidin C1", "Phloretin", "Phlorhizin", "Quercitrin")) +
  theme_minimal() +  
  theme(legend.position = "none") +
  labs(x = "Relative abundance",
       y = "Compound",
       title = "Distribution of relative abundace for some compounds in blended juice")
```

    Picking joint bandwidth of 0.563

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-24-1.png)

## Nonparametric tests

### kruskal Wallis test

``` r
kruskal_compound_B <- compound_B %>%
  drop_na(rel_abund, treatment) %>%            # Remove NA values first
  group_by(id) %>%                             # Group by 'id'
  kruskal_test(rel_abund ~ treatment) %>%      # Apply the Kruskal-Wallis test
  add_significance() %>%     
  arrange(p)     
           
kruskal_compound_B
```

    # A tibble: 6 × 8
      id             .y.           n statistic    df        p method        p.signif
      <chr>          <chr>     <int>     <dbl> <int>    <dbl> <chr>         <chr>   
    1 Procyanidin B2 rel_abund    21      18.8     3 0.000306 Kruskal-Wall… ***     
    2 Procyanidin C1 rel_abund    21      18.8     3 0.000306 Kruskal-Wall… ***     
    3 Epicatechin    rel_abund    21      18.3     3 0.000384 Kruskal-Wall… ***     
    4 Phlorhizin     rel_abund    21      17.7     3 0.000501 Kruskal-Wall… ***     
    5 Phloretin      rel_abund    21      17.1     3 0.00066  Kruskal-Wall… ***     
    6 Quercitrin     rel_abund    21      15.8     3 0.00123  Kruskal-Wall… **      

### Post-hoc test

Run Dunn test

``` r
(kruskal_compound_B_posthoc <- compound_B %>%
      drop_na(rel_abund, treatment) %>%        # Remove NA values first
      group_by(id) %>%                         # Group by 'id'
      dunn_test(rel_abund ~ treatment,
                p.adjust.method = "BH"))
```

    # A tibble: 36 × 10
       id     .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
     * <chr>  <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
     1 Epica… rel_… FJ     MEF        6     5   -2.90   3.72e-3 1.12e-2 *           
     2 Epica… rel_… FJ     MEF+SS     6     5   -4.02   5.85e-5 3.51e-4 ***         
     3 Epica… rel_… FJ     MT         6     5   -1.46   1.43e-1 2.03e-1 ns          
     4 Epica… rel_… MEF    MEF+SS     5     5   -1.07   2.85e-1 2.85e-1 ns          
     5 Epica… rel_… MEF    MT         5     5    1.38   1.69e-1 2.03e-1 ns          
     6 Epica… rel_… MEF+SS MT         5     5    2.45   1.44e-2 2.89e-2 *           
     7 Phlor… rel_… FJ     MEF        6     5   -2.10   3.55e-2 5.32e-2 ns          
     8 Phlor… rel_… FJ     MEF+SS     6     5   -2.16   3.11e-2 5.32e-2 ns          
     9 Phlor… rel_… FJ     MT         6     5    1.46   1.43e-1 1.72e-1 ns          
    10 Phlor… rel_… MEF    MEF+SS     5     5   -0.0510 9.59e-1 9.59e-1 ns          
    # ℹ 26 more rows

``` r
# Create a column called comparison

kruskal_compound_B_posthocc_1 <- kruskal_compound_B_posthoc %>%
  mutate(comparison = glue("{group1} - {group2}")) 
  
# Group data by 'id' and create a nested data frame
group_cldList <- kruskal_compound_B_posthocc_1 %>%
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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Make a table that has the maximum relative abundance for each treatment.

``` r
# Calculate maximum rel_abund for each id and treatment
max_relabund <- compound_B %>%
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
       id          treatment max_tot_relabund
       <chr>       <chr>                <dbl>
     1 Epicatechin FJ                    25.4
     2 Epicatechin MEF                   22.5
     3 Epicatechin MEF+SS                21.4
     4 Epicatechin MT                    24.8
     5 Phloretin   FJ                    17.1
     6 Phloretin   MEF                   16.3
     7 Phloretin   MEF+SS                16.3
     8 Phloretin   MT                    17.6
     9 Phlorhizin  FJ                    24.0
    10 Phlorhizin  MEF                   23.1
    # ℹ 14 more rows

Bind the groups to the maximum relative abundance.

``` r
dunn_compound_B_for_plotting <- max_relabund %>%
  left_join(group_cldList, by = c("id", "treatment"))
```

Plot

``` r
boxplot_dunn <- boxplot_compound_B +
  geom_text(data = dunn_compound_B_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = dunn_compound_B_for_plotting$Letter)) +
  labs(caption = "For each compound, treatmnets with different letters are significant different unsing Krushal Wallis test")

boxplot_dunn 
```

    Warning: Use of `dunn_compound_B_for_plotting$Letter` is discouraged.
    ℹ Use `Letter` instead.

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-30-1.png)

## Parametric tests

### ANOVA test

``` r
# Prepare your data 
clean_data <- compound_B %>%
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
       id             term         df  sumsq  meansq statistic   p.value
       <chr>          <chr>     <dbl>  <dbl>   <dbl>     <dbl>     <dbl>
     1 Phloretin      treatment     3  6.60   2.20        76.5  4.45e-10
     2 Phloretin      Residuals    17  0.489  0.0287      NA   NA       
     3 Epicatechin    treatment     3 72.6   24.2        190.   2.83e-13
     4 Epicatechin    Residuals    17  2.16   0.127       NA   NA       
     5 Phlorhizin     treatment     3 10.2    3.39        75.5  4.97e-10
     6 Phlorhizin     Residuals    17  0.764  0.0450      NA   NA       
     7 Quercitrin     treatment     3  5.04   1.68        51.6  9.43e- 9
     8 Quercitrin     Residuals    17  0.553  0.0325      NA   NA       
     9 Procyanidin B2 treatment     3 57.7   19.2        254.   2.61e-14
    10 Procyanidin B2 Residuals    17  1.29   0.0758      NA   NA       
    11 Procyanidin C1 treatment     3 37.8   12.6        257.   2.35e-14
    12 Procyanidin C1 Residuals    17  0.834  0.0491      NA   NA       

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
     1 Procyanidi… trea… MEF+SS-…          0    -3.58    -3.97     -3.20    1.78e-13
     2 Procyanidi… trea… MEF+SS-…          0    -4.33    -4.81     -3.86    1.89e-13
     3 Epicatechin trea… MEF+SS-…          0    -4.40    -5.02     -3.79    1.48e-12
     4 Procyanidi… trea… MT-MEF+…          0     3.39     2.89      3.88    2.76e-12
     5 Procyanidi… trea… MT-MEF+…          0     2.66     2.26      3.06    3.95e-12
     6 Epicatechin trea… MT-MEF+…          0     3.88     3.24      4.52    1.90e-11
     7 Epicatechin trea… MEF-FJ            0    -3.34    -3.95     -2.72    1.10e-10
     8 Procyanidi… trea… MEF-FJ            0    -2.52    -2.99     -2.05    1.56e-10
     9 Procyanidi… trea… MEF-FJ            0    -1.95    -2.33     -1.57    3.00e-10
    10 Phlorhizin  trea… MEF+SS-…          0    -1.68    -2.04     -1.31    1.56e- 9
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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Join the groups to the maximum relative abundance.

``` r
tukey_compound_B_for_plotting <- max_relabund %>%
  left_join(group_cldList_2, by = c("id", "treatment"))
```

``` r
Boxplot_tukey_compound_B <- boxplot_compound_B +
  geom_text(data = tukey_compound_B_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = Letter)) +
  labs(title = "", tag = "A")

Boxplot_tukey_compound_B
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-34-1.png)

``` r
ggsave(plot = Boxplot_tukey_compound_B,
       filename = "boxplot_compound_B_neg.svg")
```

    Saving 7 x 5 in image

# Compariosn between compounds C

## Boxplots

Filter to select the compound you want to plot.

``` r
compound_C <- feature_table_long_log2 %>%
  filter(id == "Coumaroyl quinic acid" |  id == "Feruloylquinic acid" | id == "Quercetin-3-glucoside" | id == "Sinapic acid" |id =="5-Caffeoylquinic acid" | id == "Naringin") %>%
   filter(sample_or_qc == "Sample")
```

Then plot.

``` r
boxplot_compound_C <- compound_C %>%
  ggplot(aes(x = factor(treatment, levels = c("FJ", "MT", "MEF", "MEF+SS")), 
             y = rel_abund, fill = treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_wrap(vars(id), nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
 scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5")) +
  labs(x = "Treatment",
  y = "log2(relative abundance)",
  title = "")
```

## Testing assumptions

### Normality

``` r
# testing normality by group
compound_C %>%
  drop_na(rel_abund) %>% # remove NAs
  group_by(id) %>% # test by species
  shapiro_test(rel_abund) # test for normality
```

    # A tibble: 6 × 4
      id                    variable  statistic        p
      <chr>                 <chr>         <dbl>    <dbl>
    1 5-Caffeoylquinic acid rel_abund     0.891 0.0233  
    2 Coumaroyl quinic acid rel_abund     0.830 0.00199 
    3 Feruloylquinic acid   rel_abund     0.824 0.00158 
    4 Naringin              rel_abund     0.824 0.00156 
    5 Quercetin-3-glucoside rel_abund     0.807 0.000833
    6 Sinapic acid          rel_abund     0.949 0.329   

### Constant variance

``` r
# testing constant variance

compound_C %>%
  drop_na(rel_abund) %>% # remove NAs
  levene_test(rel_abund ~ id) # test for constant variance
```

    Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    factor.

    # A tibble: 1 × 4
        df1   df2 statistic      p
      <int> <int>     <dbl>  <dbl>
    1     5   120      2.49 0.0352

``` r
compound_C %>%
  drop_na(rel_abund) %>%
  ggplot(aes(x = rel_abund, y = id, fill = id)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  # scale_y_discrete(labels = c( "Coumaroyl quinic acid", "Feruloylquinic acid", "Quercetin-3-glucoside", "Sinapic acid", "5-Caffeoylquinic acid", "Naringin")) +
  theme_minimal() +  
  theme(legend.position = "none") +
  labs(x = "Relative abundance",
       y = "Compound",
       title = "Distribution of relative abundace for some compounds in blended juice")
```

    Picking joint bandwidth of 0.236

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-40-1.png)

## Nonparametric tests

### kruskal Wallis test

``` r
kruskal_compound_C <- compound_C %>%
  drop_na(rel_abund, treatment) %>%            # Remove NA values first
  group_by(id) %>%                             # Group by 'id'
  kruskal_test(rel_abund ~ treatment) %>%      # Apply the Kruskal-Wallis test
  add_significance() %>%     
  arrange(p)     
           
kruskal_compound_C
```

    # A tibble: 6 × 8
      id                    .y.           n statistic    df        p method p.signif
      <chr>                 <chr>     <int>     <dbl> <int>    <dbl> <chr>  <chr>   
    1 5-Caffeoylquinic acid rel_abund    21      17.5     3 0.000569 Krusk… ***     
    2 Coumaroyl quinic acid rel_abund    21      17.5     3 0.000571 Krusk… ***     
    3 Naringin              rel_abund    21      16.7     3 0.00083  Krusk… ***     
    4 Quercetin-3-glucoside rel_abund    21      15.7     3 0.00128  Krusk… **      
    5 Feruloylquinic acid   rel_abund    21      15.6     3 0.00136  Krusk… **      
    6 Sinapic acid          rel_abund    21      15.1     3 0.00176  Krusk… **      

### Post-hoc test

Run Dunn test

``` r
(kruskal_compound_C_posthoc <- compound_C %>%
      drop_na(rel_abund, treatment) %>%        # Remove NA values first
      group_by(id) %>%                         # Group by 'id'
      dunn_test(rel_abund ~ treatment,
                p.adjust.method = "BH"))
```

    # A tibble: 36 × 10
       id     .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
     * <chr>  <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
     1 5-Caf… rel_… FJ     MEF        6     5    -3.17  1.54e-3 4.62e-3 **          
     2 5-Caf… rel_… FJ     MEF+SS     6     5    -3.75  1.75e-4 1.05e-3 **          
     3 5-Caf… rel_… FJ     MT         6     5    -1.46  1.43e-1 1.72e-1 ns          
     4 5-Caf… rel_… MEF    MEF+SS     5     5    -0.561 5.75e-1 5.75e-1 ns          
     5 5-Caf… rel_… MEF    MT         5     5     1.63  1.03e-1 1.54e-1 ns          
     6 5-Caf… rel_… MEF+SS MT         5     5     2.19  2.84e-2 5.68e-2 ns          
     7 Couma… rel_… FJ     MEF        6     5    -3.02  2.56e-3 7.67e-3 **          
     8 Couma… rel_… FJ     MEF+SS     6     5    -3.81  1.36e-4 8.17e-4 ***         
     9 Couma… rel_… FJ     MT         6     5    -1.37  1.72e-1 2.06e-1 ns          
    10 Couma… rel_… MEF    MEF+SS     5     5    -0.764 4.45e-1 4.45e-1 ns          
    # ℹ 26 more rows

``` r
# Create a column called comparison

kruskal_compound_C_posthocc_1 <- kruskal_compound_C_posthoc %>%
  mutate(comparison = glue("{group1} - {group2}")) 
  
# Group data by 'id' and create a nested data frame
group_cldList <- kruskal_compound_C_posthocc_1 %>%
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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Make a table that has the maximum relative abundance for each treatment.

``` r
# Calculate maximum rel_abund for each id and treatment
max_relabund <- compound_C %>%
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
       id                    treatment max_tot_relabund
       <chr>                 <chr>                <dbl>
     1 5-Caffeoylquinic acid FJ                    25.7
     2 5-Caffeoylquinic acid MEF                   25.3
     3 5-Caffeoylquinic acid MEF+SS                25.1
     4 5-Caffeoylquinic acid MT                    25.6
     5 Coumaroyl quinic acid FJ                    23.7
     6 Coumaroyl quinic acid MEF                   22.9
     7 Coumaroyl quinic acid MEF+SS                22.3
     8 Coumaroyl quinic acid MT                    23.6
     9 Feruloylquinic acid   FJ                    20.1
    10 Feruloylquinic acid   MEF                   19.6
    # ℹ 14 more rows

Bind the groups to the maximum relative abundance.

``` r
dunn_compound_C_for_plotting <- max_relabund %>%
  left_join(group_cldList, by = c("id", "treatment"))
```

Plot

``` r
boxplot_dunn <- boxplot_compound_C +
  geom_text(data = dunn_compound_C_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = dunn_compound_C_for_plotting$Letter)) +
  labs(caption = "For each compound, treatmnets with different letters are significant different unsing Krushal Wallis test")

boxplot_dunn 
```

    Warning: Use of `dunn_compound_C_for_plotting$Letter` is discouraged.
    ℹ Use `Letter` instead.

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-46-1.png)

## Parametric tests

### ANOVA test

``` r
# Prepare your data 
clean_data <- compound_C %>%
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
       id                    term         df sumsq meansq statistic  p.value
       <chr>                 <chr>     <dbl> <dbl>  <dbl>     <dbl>    <dbl>
     1 Sinapic acid          treatment     3 5.39  1.80        27.0  1.08e-6
     2 Sinapic acid          Residuals    17 1.13  0.0666      NA   NA      
     3 Coumaroyl quinic acid treatment     3 7.86  2.62        67.9  1.14e-9
     4 Coumaroyl quinic acid Residuals    17 0.657 0.0386      NA   NA      
     5 5-Caffeoylquinic acid treatment     3 1.53  0.510       37.2  1.09e-7
     6 5-Caffeoylquinic acid Residuals    17 0.233 0.0137      NA   NA      
     7 Feruloylquinic acid   treatment     3 4.92  1.64        58.4  3.67e-9
     8 Feruloylquinic acid   Residuals    17 0.477 0.0281      NA   NA      
     9 Quercetin-3-glucoside treatment     3 5.26  1.75        63.5  1.91e-9
    10 Quercetin-3-glucoside Residuals    17 0.469 0.0276      NA   NA      
    11 Naringin              treatment     3 1.85  0.617       49.4  1.32e-8
    12 Naringin              Residuals    17 0.212 0.0125      NA   NA      

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
     1 Coumaroyl … trea… MEF+SS-…          0   -1.43    -1.77     -1.10      5.33e-9
     2 Quercetin-… trea… MEF+SS-…          0   -1.18    -1.47     -0.897     7.76e-9
     3 Feruloylqu… trea… MEF+SS-…          0   -1.14    -1.43     -0.850     1.57e-8
     4 Quercetin-… trea… MT-MEF+…          0    1.16     0.856     1.45      2.14e-8
     5 Naringin    trea… MEF+SS-…          0   -0.733   -0.925    -0.540     2.68e-8
     6 Coumaroyl … trea… MT-MEF+…          0    1.32     0.962     1.67      3.75e-8
     7 Feruloylqu… trea… MT-MEF+…          0    1.11     0.806     1.41      4.58e-8
     8 Naringin    trea… MT-MEF+…          0    0.726    0.525     0.927     5.83e-8
     9 Coumaroyl … trea… MEF-FJ            0   -1.06    -1.40     -0.726     4.33e-7
    10 5-Caffeoyl… trea… MEF+SS-…          0   -0.633   -0.835    -0.432     4.37e-7
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
  
# mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Join the groups to the maximum relative abundance.

``` r
tukey_compound_C_for_plotting <- max_relabund %>%
  left_join(group_cldList_2, by = c("id", "treatment"))
```

``` r
Boxplot_tukey_compound_C <- boxplot_compound_C +
  geom_text(data = tukey_compound_C_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = Letter)) +
  labs(title = "", tag = "B")

Boxplot_tukey_compound_C
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-50-1.png)

``` r
ggsave(plot = Boxplot_tukey_compound_C,
       filename = "boxplot_compound_C_neg.svg")
```

    Saving 7 x 5 in image

# Compariosn between compounds E

## Boxplots

Filter to select the compound you want to plot.

``` r
compound_E <- feature_table_long_log2 %>%
  filter(id == "3,4-Dihydroxybenzoic acid" | id == "2,4,6-Trihydroxybenzoic acid" | id == "Glucose" | id == "2,4,6-Tryhydroxybenzaldehyde") %>%
   filter(sample_or_qc == "Sample")
```

Then plot.

``` r
boxplot_compound_E <- compound_E %>%
  ggplot(aes(x = factor(treatment, levels = c("FJ", "MT", "MEF", "MEF+SS")), 
             y = rel_abund, fill = treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3) +
  facet_wrap(vars(id), nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("FJ" = "#F28EAD", 
                                "MT" = "#C77CDB",
                                "MEF" = "#A8D8B9", 
                                "MEF+SS" = "#7FB3D5")) +
  labs(x = "Treatment",
  y = "log2(relative abundance)",
  title = "")

boxplot_compound_E
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-53-1.png)

## Testing assumptions

### Normality

``` r
# testing normality by group
compound_E %>%
  drop_na(rel_abund) %>% # remove NAs
  group_by(id) %>% # test by species
  shapiro_test(rel_abund) # test for normality
```

    # A tibble: 4 × 4
      id                           variable  statistic       p
      <chr>                        <chr>         <dbl>   <dbl>
    1 2,4,6-Trihydroxybenzoic acid rel_abund     0.848 0.00386
    2 2,4,6-Tryhydroxybenzaldehyde rel_abund     0.894 0.0266 
    3 3,4-Dihydroxybenzoic acid    rel_abund     0.878 0.0136 
    4 Glucose                      rel_abund     0.942 0.242  

### Constant variance

``` r
# testing constant variance

compound_E %>%
  drop_na(rel_abund) %>% # remove NAs
  levene_test(rel_abund ~ id) # test for constant variance
```

    Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    factor.

    # A tibble: 1 × 4
        df1   df2 statistic        p
      <int> <int>     <dbl>    <dbl>
    1     3    80      20.3 6.78e-10

``` r
compound_E %>%
  drop_na(rel_abund) %>%
  ggplot(aes(x = rel_abund, y = id, fill = id)) +
  ggridges::geom_density_ridges(alpha = 0.7) +
  # scale_y_discrete(labels = c("3,4-Dihydroxybenzoic acid", "2,4,6-Trihydroxybenzoic acid", "Glucose", "2,4,6-Trihydroxybenzaldehyde")) +
  theme_minimal() +  
  theme(legend.position = "none") +
  labs(x = "Relative abundance",
       y = "Compound",
       title = "Distribution of relative abundace for some compounds in blended juice")
```

    Picking joint bandwidth of 0.474

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-56-1.png)

## Nonparametric tests

### kruskal Wallis test

``` r
kruskal_compound_E <- compound_C %>%
  drop_na(rel_abund, treatment) %>%            # Remove NA values first
  group_by(id) %>%                             # Group by 'id'
  kruskal_test(rel_abund ~ treatment) %>%      # Apply the Kruskal-Wallis test
  add_significance() %>%     
  arrange(p)     
           
kruskal_compound_E
```

    # A tibble: 6 × 8
      id                    .y.           n statistic    df        p method p.signif
      <chr>                 <chr>     <int>     <dbl> <int>    <dbl> <chr>  <chr>   
    1 5-Caffeoylquinic acid rel_abund    21      17.5     3 0.000569 Krusk… ***     
    2 Coumaroyl quinic acid rel_abund    21      17.5     3 0.000571 Krusk… ***     
    3 Naringin              rel_abund    21      16.7     3 0.00083  Krusk… ***     
    4 Quercetin-3-glucoside rel_abund    21      15.7     3 0.00128  Krusk… **      
    5 Feruloylquinic acid   rel_abund    21      15.6     3 0.00136  Krusk… **      
    6 Sinapic acid          rel_abund    21      15.1     3 0.00176  Krusk… **      

### Post-hoc test

Run Dunn test

``` r
(kruskal_compound_E_posthoc <- compound_E %>%
      drop_na(rel_abund, treatment) %>%        # Remove NA values first
      group_by(id) %>%                         # Group by 'id'
      dunn_test(rel_abund ~ treatment,
                p.adjust.method = "BH"))
```

    # A tibble: 24 × 10
       id     .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
     * <chr>  <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
     1 2,4,6… rel_… FJ     MEF        6     5    3.15   1.64e-3 0.00491 **          
     2 2,4,6… rel_… FJ     MEF+SS     6     5    1.82   6.90e-2 0.103   ns          
     3 2,4,6… rel_… FJ     MT         6     5   -0.683  4.95e-1 0.495   ns          
     4 2,4,6… rel_… MEF    MEF+SS     5     5   -1.27   2.03e-1 0.243   ns          
     5 2,4,6… rel_… MEF    MT         5     5   -3.67   2.43e-4 0.00146 **          
     6 2,4,6… rel_… MEF+SS MT         5     5   -2.40   1.66e-2 0.0332  *           
     7 2,4,6… rel_… FJ     MEF        6     5    0.0976 9.22e-1 0.922   ns          
     8 2,4,6… rel_… FJ     MEF+SS     6     5   -2.08   3.71e-2 0.0556  ns          
     9 2,4,6… rel_… FJ     MT         6     5    2.17   2.97e-2 0.0556  ns          
    10 2,4,6… rel_… MEF    MEF+SS     5     5   -2.09   3.67e-2 0.0556  ns          
    # ℹ 14 more rows

``` r
# Create a column called comparison

kruskal_compound_E_posthocc_1 <- kruskal_compound_E_posthoc %>%
  mutate(comparison = glue("{group1} - {group2}")) 
  
# Group data by 'id' and create a nested data frame
group_cldList <- kruskal_compound_C_posthocc_1 %>%
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

 # mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Make a table that has the maximum relative abundance for each treatment.

``` r
# Calculate maximum rel_abund for each id and treatment
max_relabund <- compound_E %>%
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

    # A tibble: 16 × 3
       id                           treatment max_tot_relabund
       <chr>                        <chr>                <dbl>
     1 2,4,6-Trihydroxybenzoic acid FJ                    14.2
     2 2,4,6-Trihydroxybenzoic acid MEF                   18.7
     3 2,4,6-Trihydroxybenzoic acid MEF+SS                17.7
     4 2,4,6-Trihydroxybenzoic acid MT                    14.0
     5 2,4,6-Tryhydroxybenzaldehyde FJ                    18.2
     6 2,4,6-Tryhydroxybenzaldehyde MEF                   18.3
     7 2,4,6-Tryhydroxybenzaldehyde MEF+SS                17.6
     8 2,4,6-Tryhydroxybenzaldehyde MT                    18.6
     9 3,4-Dihydroxybenzoic acid    FJ                    18.1
    10 3,4-Dihydroxybenzoic acid    MEF                   21.4
    11 3,4-Dihydroxybenzoic acid    MEF+SS                20.6
    12 3,4-Dihydroxybenzoic acid    MT                    19.1
    13 Glucose                      FJ                    18.2
    14 Glucose                      MEF                   18.3
    15 Glucose                      MEF+SS                18.3
    16 Glucose                      MT                    18.3

Bind the groups to the maximum relative abundance.

``` r
dunn_compound_E_for_plotting <- max_relabund %>%
  left_join(group_cldList, by = c("id", "treatment"))
```

Plot

``` r
boxplot_dunn <- boxplot_compound_E +
  geom_text(data = dunn_compound_E_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = dunn_compound_E_for_plotting$Letter)) +
  labs(caption = "For each compound, treatmnets with different letters are significant different unsing Krushal Wallis test")

boxplot_dunn 
```

    Warning: Use of `dunn_compound_E_for_plotting$Letter` is discouraged.
    ℹ Use `Letter` instead.

    Warning: Removed 16 rows containing missing values (`geom_text()`).

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-62-1.png)

## Parametric tests

### ANOVA test

``` r
# Prepare your data 
clean_data <- compound_E %>%
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

    # A tibble: 8 × 7
    # Groups:   id [4]
      id                           term       df   sumsq  meansq statistic   p.value
      <chr>                        <chr>   <dbl>   <dbl>   <dbl>     <dbl>     <dbl>
    1 2,4,6-Tryhydroxybenzaldehyde treatm…     3  4.09    1.36       57.5   4.16e- 9
    2 2,4,6-Tryhydroxybenzaldehyde Residu…    17  0.404   0.0238     NA    NA       
    3 3,4-Dihydroxybenzoic acid    treatm…     3 29.9     9.96      304.    5.86e-15
    4 3,4-Dihydroxybenzoic acid    Residu…    17  0.558   0.0328     NA    NA       
    5 2,4,6-Trihydroxybenzoic acid treatm…     3 99.4    33.1       163.    9.83e-13
    6 2,4,6-Trihydroxybenzoic acid Residu…    17  3.45    0.203      NA    NA       
    7 Glucose                      treatm…     3  0.0432  0.0144      1.21  3.37e- 1
    8 Glucose                      Residu…    17  0.202   0.0119     NA    NA       

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

    # A tibble: 24 × 8
    # Groups:   id [4]
       id          term  contrast null.value estimate conf.low conf.high adj.p.value
       <chr>       <chr> <chr>         <dbl>    <dbl>    <dbl>     <dbl>       <dbl>
     1 3,4-Dihydr… trea… MEF-FJ            0     2.94    2.63       3.26    1.76e-13
     2 3,4-Dihydr… trea… MEF+SS-…          0     2.49    2.18       2.80    4.27e-13
     3 2,4,6-Trih… trea… MT-MEF            0    -5.15   -5.96      -4.34    8.24e-12
     4 2,4,6-Trih… trea… MEF-FJ            0     4.75    3.97       5.52    1.52e-11
     5 3,4-Dihydr… trea… MT-MEF            0    -1.90   -2.23      -1.58    3.31e-11
     6 2,4,6-Trih… trea… MT-MEF+…          0    -3.77   -4.58      -2.96    1.25e- 9
     7 2,4,6-Tryh… trea… MT-MEF+…          0     1.24    0.964      1.52    2.29e- 9
     8 3,4-Dihydr… trea… MT-MEF+…          0    -1.45   -1.78      -1.13    2.48e- 9
     9 2,4,6-Trih… trea… MEF+SS-…          0     3.37    2.59       4.14    3.63e- 9
    10 3,4-Dihydr… trea… MT-FJ             0     1.04    0.730      1.35    1.83e- 7
    # ℹ 14 more rows

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

#  mutate(treatment = str_replace(treatment, "FreshJuice", "Fresh Juice"))
```

Join the groups to the maximum relative abundance.

``` r
tukey_compound_E_for_plotting <- max_relabund %>%
  left_join(group_cldList_2, by = c("id", "treatment"))
```

``` r
Boxplot_tukey_compound_E <- boxplot_compound_E +
  geom_text(data = tukey_compound_E_for_plotting,
            aes(x = treatment,
                y = max_tot_relabund + 0.5,
                label = Letter)) +
  labs(caption = "Line represents the median abundance per treatment. 
       For each compound, treatmnets with different letters are significant different using Tukey's test",   
       title = "", tag = "E")

Boxplot_tukey_compound_E
```

![](Metabolomics---Neg-Data---MEF---Statistical-analyses---Part-II_files/figure-commonmark/unnamed-chunk-66-1.png)

``` r
ggsave(plot = Boxplot_tukey_compound_E,
       filename = "boxplot_compound_E_neg.svg")
```

    Saving 7 x 5 in image
