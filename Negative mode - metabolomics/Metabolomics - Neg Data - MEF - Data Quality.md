# Metabolomic first pass data analysis - Negative ionization mode - MEF
Giovana Domeneghini Mercali

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
library(rstatix) # for additional univariate functionality
library(readxl)
library(svglite)

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

## Read in data

Negative mode.

``` r
# read in metabolomics data
metab <- read_csv("Post_filtering_neg_2997_MEF.csv",
                  trim_ws = TRUE,
                  na = "0") # read in 0s to be NAs.
```

    Rows: 2997 Columns: 37
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (22): mz_rt, 40TT_10M2_neq_16, FreshJuice1_neq_24, FreshJuice3_neg_64, 4...
    dbl (15): row_ID, mz, rt, QC11_neg_78, QC4_neg_28, QC10_neg_74, QC2_neg_13, ...

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

Take a quick look at our data.

``` r
# look at first 5 rows, first 5 columns 
metab[1:5,1:5]
```

    # A tibble: 5 × 5
      row_ID mz_rt              mz    rt `40TT_10M2_neq_16`
       <dbl> <chr>           <dbl> <dbl> <chr>             
    1   1528 100.0402_0.8975  100. 0.898 12623.705         
    2    625 100.9439_0.3382  101. 0.338 NA                
    3    498 101.0222_0.3207  101. 0.321 85662.65          
    4   1055 101.0238_0.4781  101. 0.478 122471.99         
    5    345 101.0294_0.3215  101. 0.321 6636.2046         

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
dim(metab)
```

    [1] 2997   37

``` r
dim(metadata)
```

    [1] 67  6

``` r
# make metab sample columns all numeric
metab <- metab %>%
  mutate((across(.cols = 5:ncol(.), .fns = as.numeric)))
```

    Warning: There were 21 warnings in `mutate()`.
    The first warning was:
    ℹ In argument: `(across(.cols = 5:ncol(.), .fns = as.numeric))`.
    Caused by warning in `fn()`:
    ! NAs introduced by coercion
    ℹ Run `dplyr::last_dplyr_warnings()` to see the 20 remaining warnings.

``` r
glimpse(metab)
```

    Rows: 2,997
    Columns: 37
    $ row_ID             <dbl> 1528, 625, 498, 1055, 345, 1117, 1002, 630, 854, 12…
    $ mz_rt              <chr> "100.0402_0.8975", "100.9439_0.3382", "101.0222_0.3…
    $ mz                 <dbl> 100.0402, 100.9439, 101.0222, 101.0238, 101.0294, 1…
    $ rt                 <dbl> 0.8975327, 0.3381544, 0.3207379, 0.4781379, 0.32149…
    $ `40TT_10M2_neq_16` <dbl> 12623.705, NA, 85662.650, 122471.990, 6636.205, 949…
    $ FreshJuice1_neq_24 <dbl> 12023.159, NA, 88496.350, 112586.830, 6805.068, 854…
    $ FreshJuice3_neg_64 <dbl> 13652.454, NA, 96339.530, 112846.810, 7138.783, 865…
    $ `40TT_10M5_neg_40` <dbl> 11759.035, NA, 86424.410, 115492.550, 5630.096, 788…
    $ `40TT_10M4_neg_26` <dbl> 11409.089, NA, 104715.125, 119014.195, 7456.335, 94…
    $ FreshJuice4_neg_38 <dbl> 13704.267, NA, NA, 115320.030, NA, 8136.055, 6654.9…
    $ FreshJuice6_neg_45 <dbl> 11589.513, NA, 98065.090, 111706.766, 7787.599, 822…
    $ FreshJuice5_neg_73 <dbl> 13904.639, NA, 96784.640, 108390.480, 7504.919, 759…
    $ `40TT_10M6_neg_72` <dbl> 12582.365, NA, 110086.090, 123437.640, 7594.493, 94…
    $ `40TT_10M3_neg_46` <dbl> 12289.652, NA, 106494.240, 113250.050, 8251.779, 93…
    $ FreshJuice2_neg_54 <dbl> 13125.755, NA, 87831.695, 105103.160, 7000.599, 743…
    $ MEF_SS6_neg_75     <dbl> 11607.853, 13279.143, 123892.610, 110571.220, 7820.…
    $ MEF_SS4_neg_68     <dbl> 12910.261, 14138.209, 127898.520, 112346.210, 8948.…
    $ MEF_SS1_neg_56     <dbl> 10724.819, 16645.283, 131427.330, 103126.720, 10169…
    $ MEF_SS2_neg_34     <dbl> 9979.047, 15065.035, 130573.820, 109426.490, 8826.8…
    $ MEF_SS3_neg_52     <dbl> 11184.081, 13860.484, 128510.580, 115518.510, 8788.…
    $ MEF1_neg_41        <dbl> 10983.453, 96104.220, 108690.195, 108788.445, 7488.…
    $ MEF2_neg_60        <dbl> 11319.180, 43299.195, 112216.840, 109649.630, 7470.…
    $ QC11_neg_78        <dbl> 13470.114, 5929.262, 120178.180, 125090.984, 8608.3…
    $ QC4_neg_28         <dbl> 13730.859, 5456.804, 108522.360, 118438.914, 7660.1…
    $ QC10_neg_74        <dbl> 12069.027, 4993.207, 109453.290, 118015.520, 7795.1…
    $ QC2_neg_13         <dbl> 13325.600, 6165.458, 119287.766, 125860.680, 9659.3…
    $ MEF4_neg_29        <dbl> 11358.533, 19994.375, 132241.270, 93996.280, 8747.0…
    $ MEF5_neg_50        <dbl> 11670.782, 21388.893, 116037.100, 111779.020, 8277.…
    $ MEF3_neg_44        <dbl> 13886.385, 28309.137, 117381.240, 109532.600, 8700.…
    $ QC9_neg_67         <dbl> 11914.959, 6696.042, 92050.380, 113209.550, 6111.62…
    $ QC5_neg_35         <dbl> 11789.217, 7697.496, 101195.160, 112213.750, 6851.1…
    $ QC6_neg_43         <dbl> 11772.166, 4621.732, 176979.610, 116165.266, 11990.…
    $ QC3_neq_20         <dbl> 11308.223, 3726.545, 119111.340, 121476.016, 9253.6…
    $ QC7_neg_51         <dbl> 11150.521, 3754.003, 109563.060, 115072.195, 7161.4…
    $ QC8_neg_59         <dbl> 11563.439, 6900.274, 114135.940, 119587.280, 8842.6…
    $ QC12_neg_79        <dbl> 12009.428, 5218.626, 112950.984, 112701.120, 7004.2…
    $ QC1_neg_12         <dbl> 11507.536, 7507.896, 119697.086, 129105.930, 8513.6…

## Wrangle sample names

Here, the samples are in columns and the features are in rows. Samples
are coded so that the first number is the treatment code, and the last
code is the run order. We are going to transpose the data so that
samples are in rows and features are in columns, and we will also import
the metadata about the samples.

``` r
metab_t <- metab %>%
  select(-row_ID, -mz, -rt) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "mz_rt")

# make the first row the column names
colnames(metab_t) <- metab_t[1,]

# then remove the first row and rename the column that should be called sample_name
metab_t <- metab_t[-1,] %>%
  rename(sample_name = mz_rt) %>%
  mutate((across(.cols = 2:ncol(.), .fns = as.numeric)))

metab_t[1:10, 1:10]
```

              sample_name 100.0402_0.8975 100.9439_0.3382 101.0222_0.3207
    2    40TT_10M2_neq_16        12623.70              NA        85662.65
    3  FreshJuice1_neq_24        12023.16              NA        88496.35
    4  FreshJuice3_neg_64        13652.45              NA        96339.53
    5    40TT_10M5_neg_40        11759.03              NA        86424.41
    6    40TT_10M4_neg_26        11409.09              NA       104715.12
    7  FreshJuice4_neg_38        13704.27              NA              NA
    8  FreshJuice6_neg_45        11589.51              NA        98065.09
    9  FreshJuice5_neg_73        13904.64              NA        96784.64
    10   40TT_10M6_neg_72        12582.36              NA       110086.09
    11   40TT_10M3_neg_46        12289.65              NA       106494.24
       101.0238_0.4781 101.0294_0.3215 101.031_0.4766 103.003_0.4435
    2         122472.0        6636.205       9494.004       4520.221
    3         112586.8        6805.068       8547.527       5503.039
    4         112846.8        7138.783       8658.020       6064.336
    5         115492.6        5630.096       7888.324       5544.878
    6         119014.2        7456.335       9419.192       5130.791
    7         115320.0              NA       8136.055       6654.931
    8         111706.8        7787.599       8229.951       6015.914
    9         108390.5        7504.919       7593.521       5728.873
    10        123437.6        7594.493       9404.909       5407.939
    11        113250.1        8251.779       9371.131       5644.169
       103.0384_0.3623 103.0396_0.3544
    2         14905.29        19076.14
    3         16517.44        18924.65
    4         17539.24        18406.06
    5         14776.22        23682.22
    6         19719.17        21650.73
    7         17290.20        18125.31
    8         15973.13        17754.90
    9         16051.92        19278.06
    10        20082.85        23637.32
    11        22339.99        22790.52

Add in the metadata and make a new column that will indicate whether a
sample is a “sample” or a “QC”. The metadata we have for the samples we
are no longer using (process blanks etc) are removed.

``` r
metab_plus <- left_join(metab_t, metadata, by = "sample_name") %>% # keeps only those samples in metab_t
  mutate(sample_or_qc = if_else(str_detect(sample_name, "QC"), true = "QC", false = "Sample")) %>%
  dplyr::select(sample_name, short_sample_name, treatment, run_order, soluble_solids, mass, sample_or_qc, everything()) %>% # and move metadata to the front
mutate((across(.cols = 5:6, .fns = as.numeric)))
```

Go from wide to long data.

``` r
metab_plus_long <- metab_plus %>%
  pivot_longer(cols = -c(sample_name, short_sample_name, treatment, run_order, sample_or_qc, soluble_solids, mass),  # remove metadata
               names_to = "mz_rt",
               values_to = "rel_abund")

glimpse(metab_plus_long)
```

    Rows: 98,901
    Columns: 9
    $ sample_name       <chr> "40TT_10M2_neq_16", "40TT_10M2_neq_16", "40TT_10M2_n…
    $ short_sample_name <chr> "MT-2", "MT-2", "MT-2", "MT-2", "MT-2", "MT-2", "MT-…
    $ treatment         <chr> "MT", "MT", "MT", "MT", "MT", "MT", "MT", "MT", "MT"…
    $ run_order         <dbl> 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, …
    $ soluble_solids    <dbl> 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3…
    $ mass              <dbl> 0.4284, 0.4284, 0.4284, 0.4284, 0.4284, 0.4284, 0.42…
    $ sample_or_qc      <chr> "Sample", "Sample", "Sample", "Sample", "Sample", "S…
    $ mz_rt             <chr> "100.0402_0.8975", "100.9439_0.3382", "101.0222_0.32…
    $ rel_abund         <dbl> 12623.705, NA, 85662.650, 122471.990, 6636.205, 9494…

Also add separate columns for mz and rt, and making both numeric.

``` r
metab_plus_long <- metab_plus_long %>%
  separate_wider_delim(cols = mz_rt,
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE) %>%
  mutate(across(.cols = 8:9, .fns = as.numeric))  %>% # convert mz and rt to numeric
  drop_na(mz, rt) # removing mz and rt values that are NA; maybe not actually necessary

# how did that go?
head(metab_plus_long)
```

    # A tibble: 6 × 11
      sample_name      short_sample_name treatment run_order soluble_solids  mass
      <chr>            <chr>             <chr>         <dbl>          <dbl> <dbl>
    1 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    2 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    3 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    4 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    5 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    6 40TT_10M2_neq_16 MT-2              MT               16           10.3 0.428
    # ℹ 5 more variables: sample_or_qc <chr>, mz <dbl>, rt <dbl>, mz_rt <chr>,
    #   rel_abund <dbl>

## Correction - soluble solids and mass

``` r
metab_corrected_long <- metab_plus_long |>
 mutate(rel_abund = if_else(!is.na(rel_abund), true = rel_abund*40/(soluble_solids*mass), false = rel_abund)) |>
  select(-soluble_solids, -mass, -sample_name) |>
  rename("sample_name" = short_sample_name)

head(metab_corrected_long)
```

    # A tibble: 6 × 8
      sample_name treatment run_order sample_or_qc    mz    rt mz_rt       rel_abund
      <chr>       <chr>         <dbl> <chr>        <dbl> <dbl> <chr>           <dbl>
    1 MT-2        MT               16 Sample        100. 0.898 100.0402_0…   114435.
    2 MT-2        MT               16 Sample        101. 0.338 100.9439_0…       NA 
    3 MT-2        MT               16 Sample        101. 0.321 101.0222_0…   776542.
    4 MT-2        MT               16 Sample        101. 0.478 101.0238_0…  1110223.
    5 MT-2        MT               16 Sample        101. 0.322 101.0294_0…    60158.
    6 MT-2        MT               16 Sample        101. 0.477 101.031_0.…    86064.

## Data summaries

What mass range do I have?

``` r
range(metab_corrected_long$mz)
```

    [1]  100.0402 1509.3982

What retention time range do I have?

``` r
range(metab_corrected_long$rt)
```

    [1] 0.3003 7.4500

How many samples are in each of my meta-data groups?

``` r
# make wide data to make some calculations easier
metab_corrected_wide <- metab_corrected_long %>%
  dplyr::select(-mz, -rt) %>%
  pivot_wider(names_from = mz_rt,
              values_from = rel_abund)

# by sample vs QC
metab_corrected_wide %>%
  count(sample_or_qc)
```

    # A tibble: 2 × 2
      sample_or_qc     n
      <chr>        <int>
    1 QC              12
    2 Sample          21

What does my data coverage across mz and rt look like?

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot (all features)",
       subtitle = "C18 reversed phase, negative ionization mode")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-13-1.png)

Distribution of masses

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>%
  ggplot(aes(x = mz)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  labs(x = "Mass to charge ratio (m/z)",
       y = "Number of features",
       title = "Distribution of features by mass")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-14-1.png)

Distribution of retention times

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>%
  ggplot(aes(x = rt)) +
  geom_density() +
  theme_minimal() +
  labs(x = "Retention time",
       y = "Number of features",
       title = "Distribution of features by retention time")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-15-1.png)

## Missing data

### Surveying missingness

How many missing values are there for each feature? In this dataset,
missing values are coded as zero.

``` r
# all data including QCs
# how many missing values are there for each feature (row)
na_by_feature <- rowSums(is.na(metab)) %>%
  as.data.frame() %>%
  rename(missing_values = 1)


na_by_feature %>%
  ggplot(aes(x = missing_values)) +
  geom_histogram(bins = 67) + # since 67 samples
  theme_minimal() + 
  labs(title = "Number of missing values for each feature",
       x = "Number of missing values",
       y = "How many features have that many missing values")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-16-1.png)

How many features have no missing values?

``` r
na_by_feature %>%
  count(missing_values == 0)
```

      missing_values == 0    n
    1               FALSE  389
    2                TRUE 2608

How many missing values are there for each sample?

``` r
# all data including QCs
# how many missing values are there for each feature (row)
na_by_sample <- colSums(is.na(metab)) %>%
  as.data.frame() %>%
  rename(missing_values = 1) %>%
  rownames_to_column(var = "feature") %>%
  filter(!feature %in% c("mz_rt", "rt", "mz", "row_ID"))

na_by_sample %>%
  ggplot(aes(x = missing_values)) +
  geom_histogram(bins = 100) + # there are 3473 features
  theme_minimal() + 
  labs(title = "Number of missing values for each sample",
       x = "Number of missing values",
       y = "How many samples have that many missing values")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-18-1.png)

Which features have a lot of missing values?

``` r
contains_NAs_feature <- metab_corrected_long %>%
  group_by(mz_rt) %>%
  count(is.na(rel_abund)) %>%
  filter(`is.na(rel_abund)` == TRUE) %>% 
  arrange(desc(n))

head(contains_NAs_feature)
```

    # A tibble: 6 × 3
    # Groups:   mz_rt [6]
      mz_rt           `is.na(rel_abund)`     n
      <chr>           <lgl>              <int>
    1 419.0805_4.825  TRUE                  21
    2 636.5316_6.7983 TRUE                  21
    3 708.4612_6.7979 TRUE                  21
    4 419.067_4.822   TRUE                  19
    5 421.0648_4.8238 TRUE                  17
    6 227.9899_3.5965 TRUE                  16

Which samples have a lot of missing values?

``` r
contains_NAs_sample <- metab_corrected_long %>%
  group_by(sample_name) %>%
  count(is.na(rel_abund)) %>%
  filter(`is.na(rel_abund)` == TRUE) %>%
  arrange(desc(n))

head(contains_NAs_sample)
```

    # A tibble: 6 × 3
    # Groups:   sample_name [6]
      sample_name `is.na(rel_abund)`     n
      <chr>       <lgl>              <int>
    1 FJ-4        TRUE                 195
    2 FJ-1        TRUE                 194
    3 FJ-5        TRUE                 192
    4 FJ-2        TRUE                 188
    5 FJ-6        TRUE                 183
    6 FJ-3        TRUE                 180

Are there any missing values in the QCs? (There shouldn’t be.)

``` r
metab_QC <- metab %>%
  dplyr::select(contains("QC"))

na_by_sample <- colSums(is.na(metab_QC)) %>%
  as.data.frame() %>%
  rename(missing_values = 1) %>%
  rownames_to_column(var = "feature") %>%
  filter(!feature == "mz_rt")

sum(na_by_sample$missing_values) # nope
```

    [1] 0

### Imputing missing values

This is an optional step but some downstream analyses don’t handle
missingness well. Here we are imputing missing data with half the lowest
value observed for that feature.

``` r
# grab only the feature data and leave metadata
metab_corrected_wide_imputed <- metab_corrected_wide %>%
  dplyr::select(-c(1:4)) # the metadata columns

metab_corrected_wide_imputed[] <- lapply(metab_corrected_wide_imputed,
                                 function(x) ifelse(is.na(x), min(x, na.rm = TRUE)/2, x))

metab_corrected_wide_imputed_notame <- metab_corrected_wide_imputed

# bind back the metadata
metab_corrected_wide_imputed <- bind_cols(metab_corrected_wide[,1:4], metab_corrected_wide_imputed)
```

For notame later:

``` r
data_notame <- bind_cols(metab_corrected_wide[,1], metab_corrected_wide_imputed_notame) %>%
 as.data.frame()
```

Did imputing work?

``` r
# count missing values
metab_corrected_wide_imputed %>%
  dplyr::select(-c(1:4)) %>% # where the metadata is
  is.na() %>%
  sum()
```

    [1] 0

Create long imputed dataset.

``` r
metab_corrected_long_imputed <- metab_corrected_wide_imputed %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund")

head(metab_corrected_long_imputed)
```

    # A tibble: 6 × 6
      sample_name treatment run_order sample_or_qc mz_rt           rel_abund
      <chr>       <chr>         <dbl> <chr>        <chr>               <dbl>
    1 MT-2        MT               16 Sample       100.0402_0.8975   114435.
    2 MT-2        MT               16 Sample       100.9439_0.3382    17707.
    3 MT-2        MT               16 Sample       101.0222_0.3207   776542.
    4 MT-2        MT               16 Sample       101.0238_0.4781  1110223.
    5 MT-2        MT               16 Sample       101.0294_0.3215    60158.
    6 MT-2        MT               16 Sample       101.031_0.4766     86064.

Let’s also make separate mz and rt columns.

``` r
metab_corrected_long_imputed <- metab_corrected_long_imputed %>%
  separate_wider_delim(cols = mz_rt,
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE)

metab_corrected_long_imputed$mz <- as.numeric(metab_corrected_long_imputed$mz)
metab_corrected_long_imputed$rt <- as.numeric(metab_corrected_long_imputed$rt)
```

## Feature clustering with `notame`

We want to cluster features that likely come from the same metabolite
together, and we can do this using the package `notame`. You can learn
more
[here](http://127.0.0.1:24885/library/notame/doc/feature_clustering.html).

``` r
browseVignettes("notame")
```

Let’s make a m/z by retention time plot again before we start.

``` r
(before_notame <- metab_corrected_long_imputed %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot before notame, 2997 features",
       subtitle = "C18 reverse phase, negative ionization mode"))
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-28-1.png)

### Wrangling data

Transpose the wide data for notame and wrangle to the right format.
Below is info from the documentation:

- Data should be a data frame containing the abundances of features in
  each sample, one row per sample, each feature in a separate column
- Features should be a data frame containing information about the
  features, namely feature name (should be the same as the column name
  in data), mass and retention time

``` r
# change to results to numeric
data_notame <- data_notame %>%
  mutate(across(-sample_name, as.numeric))

tibble(data_notame)
```

    # A tibble: 33 × 2,998
       sample_name `100.0402_0.8975` `100.9439_0.3382` `101.0222_0.3207`
       <chr>                   <dbl>             <dbl>             <dbl>
     1 MT-2                  114435.            17707.           776542.
     2 FJ-1                  108284.            17707.           797020.
     3 FJ-3                  120828.            17707.           852630.
     4 MT-5                  107273.            17707.           788416.
     5 MT-4                  103128.            17707.           946535.
     6 FJ-4                  121842.            17707.           388271.
     7 FJ-6                  105976.            17707.           896716.
     8 FJ-5                  126520.            17707.           880653.
     9 MT-6                  115773.            17707.          1012929.
    10 MT-3                  111484.            17707.           966052.
    # ℹ 23 more rows
    # ℹ 2,994 more variables: `101.0238_0.4781` <dbl>, `101.0294_0.3215` <dbl>,
    #   `101.031_0.4766` <dbl>, `103.003_0.4435` <dbl>, `103.0384_0.3623` <dbl>,
    #   `103.0396_0.3544` <dbl>, `103.0397_0.6197` <dbl>, `103.0399_0.8555` <dbl>,
    #   `103.0399_1.01` <dbl>, `107.0348_0.9384` <dbl>, `109.0295_1.8499` <dbl>,
    #   `111.0082_0.5501` <dbl>, `111.0087_0.3485` <dbl>, `111.0087_1.3341` <dbl>,
    #   `111.0087_2.2117` <dbl>, `111.0088_1.9584` <dbl>, …

Create df with features.

``` r
features <- metab_corrected_long_imputed %>%
  dplyr::select(mz_rt, mz, rt) %>%
  mutate(across(c(mz, rt), as.numeric)) %>%
  as.data.frame() %>%
  distinct()

glimpse(features)
```

    Rows: 2,997
    Columns: 3
    $ mz_rt <chr> "100.0402_0.8975", "100.9439_0.3382", "101.0222_0.3207", "101.02…
    $ mz    <dbl> 100.0402, 100.9439, 101.0222, 101.0238, 101.0294, 101.0310, 103.…
    $ rt    <dbl> 0.8975, 0.3382, 0.3207, 0.4781, 0.3215, 0.4766, 0.4435, 0.3623, …

``` r
class(features)
```

    [1] "data.frame"

### Find connections

Set `cache = TRUE` for this chunk since its a bit slow especially if you
have a lot of features (this step took ~10 min).

``` r
connection <- find_connections(data = data_notame,
                               features = features,
                               corr_thresh = 0.9,
                               rt_window = 1/60,
                               name_col = "mz_rt",
                               mz_col = "mz",
                               rt_col = "rt")
```

    Warning: executing %dopar% sequentially: no parallel backend registered

    [1] 100
    [1] 200
    [1] 300
    [1] 400
    [1] 500
    [1] 600
    [1] 700
    [1] 800
    [1] 900
    [1] 1000
    [1] 1100
    [1] 1200
    [1] 1300
    [1] 1400
    [1] 1500
    [1] 1600
    [1] 1700
    [1] 1800
    [1] 1900
    [1] 2000
    [1] 2100
    [1] 2200
    [1] 2300
    [1] 2400
    [1] 2500
    [1] 2600
    [1] 2700
    [1] 2800
    [1] 2900

``` r
head(connection)
```

                    x               y       cor rt_diff  mz_diff
    1 101.0222_0.3207 101.0294_0.3215 0.9438452  0.0008   0.0072
    2 109.0295_1.8499 153.0194_1.8489 0.9663899 -0.0010  43.9899
    3 109.0295_1.8499 153.0281_1.8495 0.9587234 -0.0004  43.9986
    4 109.0295_1.8499 377.1084_1.8398 0.9421770 -0.0101 268.0789
    5 111.0087_0.3485 111.0152_0.3494 0.9186789  0.0009   0.0065
    6 111.0087_0.3485 174.0215_0.3497 0.9011520  0.0012  63.0128

### Clustering

Now that we have found all of the features that are connected based on
the parameters we have set, we now need to find clusters.

``` r
clusters <- find_clusters(connections = connection, 
                          d_thresh = 0.8)
```

    187 components found

    Warning: executing %dopar% sequentially: no parallel backend registered

    Component 100 / 187 
    157 components found

    Component 100 / 157 
    72 components found

    27 components found

    9 components found

    7 components found

    5 components found

    5 components found

    1 components found

Assign a cluster ID to each feature to keep, and the feature that is
picked is the one with the highest median peak intensity across the
samples.

``` r
# assign a cluster ID to all features
# clusters are named after feature with highest median peak height
features_clustered <- assign_cluster_id(data_notame, 
                                        clusters, 
                                        features, 
                                        name_col = "mz_rt")
```

Export out a list of your clusters this way you can use this later
during metabolite ID.

``` r
# export clustered feature list this way you have it
write_csv(features_clustered,
          "features_clustered_neg_MEF.csv")
```

Pull data out from the clusters and see how many features we
removed/have now.

``` r
# lets see how many features are removed when we only keep one feature per cluster
pulled <- pull_clusters(data_notame, features_clustered, name_col = "mz_rt")

cluster_data <- pulled$cdata
cluster_features <- pulled$cfeatures

# how many features did we originally have after filtering?
nrow(metab)
```

    [1] 2997

``` r
# how many features got removed during clustering?
nrow(metab) - nrow(cluster_features)
```

    [1] 712

``` r
# what percentage of the original features were removed?
((nrow(metab) - nrow(cluster_features))/nrow(metab)) * 100
```

    [1] 23.75709

Reduce our dataset to include only our new clusters. `cluster_data`
contains only the retained clusters, while `cluster_features` tells you
also which features are a part of each cluster.

``` r
# combined metadata_plus with cluster_features

metab_imputed_clustered_wide <- left_join(metab_corrected_wide_imputed[,1:4], cluster_data,
                                          by = "sample_name") 

dim(metab_imputed_clustered_wide) # we have 2474 features since 4 metadata columns
```

    [1]   33 2289

``` r
# make a long/tidy df
metab_imputed_clustered_long <- metab_imputed_clustered_wide %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund") %>%
  separate_wider_delim(cols = mz_rt, # make separate columns for mz and rt too
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE) %>%
  mutate(across(.cols = c("mz", "rt"), .fns = as.numeric)) # make mz and rt numeric
```

Write out the final clustered dataset.

``` r
write_csv(metab_imputed_clustered_long,
  "features_clustered_long__neg_MEF.csv")
```

Let’s look at a m/z by retention time plot again after clustering.

``` r
(after_notame <- metab_imputed_clustered_long %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot after notame, 2289 features",
       subtitle = "C18 reverse phase, negative ionization mode"))
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-39-1.png)

``` r
before_notame / after_notame
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-40-1.png)

## Assessing data quality

Let’s make sure that our data is of good quality.

### Untransformed data

First we are going to convert the type of some of the columns to match
what we want (e.g., run order converted to numeric, treatment to
factor).

``` r
tibble(metab_imputed_clustered_long)
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
metab_imputed_clustered_long$run_order <- as.numeric(metab_imputed_clustered_long$run_order)

# make treatment and sample_or_qc a factor (i.e., categorical)
metab_imputed_clustered_long$treatment <- as.factor(metab_imputed_clustered_long$treatment)
metab_imputed_clustered_long$sample_or_qc <- as.factor(metab_imputed_clustered_long$sample_or_qc)

# did it work?
tibble(metab_imputed_clustered_long)
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

Let’s make a boxplot to see how the metabolite abundance looks across
different samples.

``` r
metab_imputed_clustered_long %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is unscaled",
       x = "Sample",
       y = "Relative abundance")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-42-1.png)

Can’t really see anything because data is skewed.

### Transformed data

#### Log2 transformed

We will log2 transform our data.

``` r
metab_imputed_clustered_long_log2 <- metab_imputed_clustered_long %>%
  mutate(rel_abund = log2(rel_abund))
```

And then plot.

``` r
boxplot_log2_transformed <- metab_imputed_clustered_long_log2 %>%
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

boxplot_log2_transformed 
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-44-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
boxplot_log2_transformed_run_order <- metab_imputed_clustered_long_log2 %>%
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

boxplot_log2_transformed_run_order
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-45-1.png)

#### Log10 transformed

We will log10 transform our data.

``` r
metab_imputed_clustered_long_log10 <- metab_imputed_clustered_long %>%
  mutate(rel_abund = log10(rel_abund))
```

We can look at this data where we group by species.

``` r
metab_imputed_clustered_long_log10 %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance", x = "Sample")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-47-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
metab_imputed_clustered_long_log10 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (+) feature abundances by sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance", x = "Sample")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-48-1.png)

#### Autoscaled

Scales to unit variance, where the mean for each feature is 0 with a
standard deviation of 1. This works well when all metabolites are
considered to be equivalently important though measurement errors can be
inflated. We’ve never actually scaled data this way for a project but
whatever here it is.

``` r
metab_imputed_clustered_wide <- metab_imputed_clustered_long %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")

# autoscaling on the now zero-centered matrix
autoscaled <- 
    bind_cols(metab_imputed_clustered_wide[1:4], # metadata
            lapply(metab_imputed_clustered_wide[5:ncol(metab_imputed_clustered_wide)], # metab data
                   scale)) # scale to mean 0 sd 1

autoscaled[1:10,1:10]
```

    # A tibble: 10 × 10
       sample_name treatment run_order sample_or_qc `100.0402_0.8975`[,1]
       <chr>       <fct>         <dbl> <fct>                        <dbl>
     1 MT-2        MT               16 Sample                     -0.101 
     2 FJ-1        FJ               24 Sample                     -0.878 
     3 FJ-3        FJ               64 Sample                      0.706 
     4 MT-5        MT               40 Sample                     -1.01  
     5 MT-4        MT               26 Sample                     -1.53  
     6 FJ-4        FJ               38 Sample                      0.834 
     7 FJ-6        FJ               45 Sample                     -1.17  
     8 FJ-5        FJ               73 Sample                      1.43  
     9 MT-6        MT               72 Sample                      0.0677
    10 MT-3        MT               46 Sample                     -0.474 
    # ℹ 5 more variables: `100.9439_0.3382` <dbl[,1]>, `101.0222_0.3207` <dbl[,1]>,
    #   `101.0238_0.4781` <dbl[,1]>, `101.031_0.4766` <dbl[,1]>,
    #   `103.003_0.4435` <dbl[,1]>

``` r
autoscaled_long <- autoscaled %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund") %>%
  mutate(treatment = as.factor(treatment))

autoscaled_long %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) Feature Abundances by Sample",
       subtitle = "Data is autoscaled",
       y = "Relative abundance")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-50-1.png)

This is weird we won’t use this.

#### Pareto scaled

Pareto scaling scales but keeps the fidelity of the original differences
in absolute measurement value more than autoscaling. Often data is
Pareto scaled after log transformation

``` r
metab_imputed_clustered_wide_log2 <- metab_imputed_clustered_long_log2 %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")

metab_imputed_clustered_wide_log2_metabs <- 
  metab_imputed_clustered_wide_log2[,5:ncol(metab_imputed_clustered_wide_log2)]

pareto_scaled <- 
  IMIFA::pareto_scale(metab_imputed_clustered_wide_log2_metabs, center = FALSE)

pareto_scaled <- bind_cols(metab_imputed_clustered_wide_log2[,1:4], pareto_scaled)
```

``` r
pareto_scaled_long <- pareto_scaled %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund")

pareto_scaled_long %>%
  # mutate(short_sample_name = fct_reorder(short_sample_name, treatment)) %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is Pareto scaled",
       y = "Relative abundance", x = "Sample")
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-51-1.png)

I think pareto scaling is making everything look super the same. I am
going to use log2 transformed data for the rest of this analysis.

Let’s write out that data as our final feature table used for the rest
of the analysis.

``` r
metab_imputed_clustered_wide_log2 <- metab_imputed_clustered_long_log2 %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = mz_rt,
              values_from = rel_abund)

write_csv(metab_imputed_clustered_wide_log2,
          file = "features_clustered_wide_log2_neg_MEF.csv")
```

## PCAs all data

### With QCs

``` r
pca_qc <- prcomp(metab_imputed_clustered_wide_log2[,-c(1:4)], # remove metadata
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
    Standard deviation     1.37087 1.35399 1.33812 1.2932 4.668e-14
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
    Standard deviation     1.338117 1.293193 4.668325e-14
    Proportion of Variance 0.002040 0.001900 0.000000e+00
    Cumulative Proportion  0.998100 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-55-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-56-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_qc <- as.data.frame(pca_qc$x)

# bind meta-data
scores_qc <- bind_cols(metab_imputed_clustered_wide_log2[,1:4], # first 4 columns
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
  geom_point(shape = 21, color = "black") +
  theme_minimal() +
  labs(x = glue("PC1: {PC1_percent_qc}%"), 
       y = glue("PC2: {PC2_percent_qc}%"), 
       title = "PCA scores plot colored by treatment, negative mode"))
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-58-1.png)

Then make your scores plot ineractive so you can see who is who.

``` r
ggplotly(scores_qc_plot, tooltip = "text")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-60-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

### Without QCs

``` r
metab_imputed_clustered_wide_log2_noqc <- metab_imputed_clustered_wide_log2 %>%
  filter(sample_or_qc == "Sample")


pca_noqc <- prcomp(metab_imputed_clustered_wide_log2_noqc[,-c(1:4)], # remove metadata
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
                              PC15    PC16    PC17    PC18    PC19    PC20     PC21
    Standard deviation     2.16276 2.12959 2.12624 2.06323 1.98795 1.86030 4.56e-14
    Proportion of Variance 0.00368 0.00356 0.00355 0.00335 0.00311 0.00272 0.00e+00
    Cumulative Proportion  0.98371 0.98728 0.99083 0.99417 0.99728 1.00000 1.00e+00

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
    Standard deviation     1.98795 1.860302 4.560158e-14
    Proportion of Variance 0.00311 0.002720 0.000000e+00
    Cumulative Proportion  0.99728 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-63-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-64-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_noqc <- as.data.frame(pca_noqc$x)

# bind meta-data
scores_noqc <- bind_cols(metab_imputed_clustered_wide_log2_noqc[,1:4], # metadata
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
  geom_point(shape = 21, color = "black") +
  theme_minimal() +
  labs(x = glue("PC1: {PC1_percent_noqc}%"), 
       y = glue("PC2: {PC2_percent_noqc}%"), 
       title = "PCA scores plot colored by treatment, negative mode"))
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-66-1.png)

Then make your scores plot ineractive so you can see who is who.

``` r
ggplotly(scores_noqc_plot, tooltip = "text")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_noqc)
```

![](Metabolomics---Neg-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-68-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.
