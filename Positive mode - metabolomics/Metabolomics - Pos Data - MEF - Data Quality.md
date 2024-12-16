# Metabolomic first pass data analysis - Positive ionization mode - MEF
Giovana Domeneghini Mercali

<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/htmlwidgets-1.6.2/htmlwidgets.js"></script>
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/plotly-binding-4.10.2/plotly.js"></script>
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/setprototypeof-0.1/setprototypeof.js"></script>
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/typedarray-0.1/typedarray.min.js"></script>
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/jquery-3.5.1/jquery.min.js"></script>
<link href="Metabolomics - Pos Data - MEF - Data Quality_files/libs/crosstalk-1.2.0/css/crosstalk.min.css" rel="stylesheet" />
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/crosstalk-1.2.0/js/crosstalk.min.js"></script>
<link href="Metabolomics - Pos Data - MEF - Data Quality_files/libs/plotly-htmlwidgets-css-2.11.1/plotly-htmlwidgets.css" rel="stylesheet" />
<script src="Metabolomics - Pos Data - MEF - Data Quality_files/libs/plotly-main-2.11.1/plotly-latest.min.js"></script>


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

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

## Read in data

Positive mode.

``` r
# read in metabolomics data
metab <- read_csv("Post_filtering_pos_3289_MEF.csv",
                  trim_ws = TRUE,
                  na = "0") # read in 0s to be NAs.
```

    Rows: 3289 Columns: 37
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (22): mz_rt, FreshJuice5_pos_73, FreshJuice3_pos_64, 40TT_10M4_pos_26, F...
    dbl (15): row_ID, mz, rt, QC6_pos_43, QC9_pos_67, QC8_pos_59, QC10_pos_74, Q...

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

Take a quick look at our data.

``` r
# look at first 5 rows, first 5 columns 
metab[1:5,1:5]
```

    # A tibble: 5 × 5
      row_ID mz_rt              mz    rt FreshJuice5_pos_73
       <dbl> <chr>           <dbl> <dbl> <chr>             
    1    906 100.0394_0.4966  100. 0.497 NA                
    2    246 100.0748_0.3036  100. 0.304 11359.0625        
    3   2770 100.112_1.8034   100. 1.80  9083.935          
    4   2472 100.1121_1.5856  100. 1.59  NA                
    5   2629 101.0232_1.6779  101. 1.68  5578.2925         

``` r
# look at first 5 rows, all columns 
metadata[1:5,]
```

    # A tibble: 5 × 6
      sample_name        short_sample_name treatment run_order soluble_solids  mass
      <chr>              <chr>             <chr>         <dbl>          <dbl> <dbl>
    1 FreshJuice2_pos_54 FJ-2              FJ               54           10.3 0.428
    2 40TT_10M6_pos_72   MT-6              MT               72           10.2 0.426
    3 40TT4_pos_49       40TT-4            40TT             49            9.9 0.425
    4 40TT1_pos_53       40TT-1            40TT             53           10.1 0.434
    5 40TT6_pos_69       40TT-6            40TT             69           10   0.429

``` r
# check dimensions
dim(metab)
```

    [1] 3289   37

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

    Rows: 3,289
    Columns: 37
    $ row_ID             <dbl> 906, 246, 2770, 2472, 2629, 750, 2370, 1895, 1299, …
    $ mz_rt              <chr> "100.0394_0.4966", "100.0748_0.3036", "100.112_1.80…
    $ mz                 <dbl> 100.0394, 100.0748, 100.1120, 100.1121, 101.0232, 1…
    $ rt                 <dbl> 0.4966192, 0.3035967, 1.8034459, 1.5856134, 1.67791…
    $ FreshJuice5_pos_73 <dbl> NA, 11359.062, 9083.935, NA, 5578.292, 53138.734, 5…
    $ FreshJuice3_pos_64 <dbl> NA, 6991.805, 10544.096, NA, 7092.719, 61444.938, 5…
    $ `40TT_10M4_pos_26` <dbl> NA, 9396.625, 7078.440, NA, 14991.171, 64184.207, 1…
    $ FreshJuice4_pos_38 <dbl> NA, 7958.649, 9853.188, 1311.417, 5236.963, 49287.3…
    $ FreshJuice6_pos_45 <dbl> NA, 9622.346, 9298.826, NA, 4703.250, 47034.375, 47…
    $ FreshJuice1_pos_24 <dbl> NA, 7768.952, 8173.136, NA, 3129.431, 51654.850, 47…
    $ `40TT_10M3_pos_46` <dbl> NA, 7989.313, 7782.835, NA, 19899.783, 67361.870, 1…
    $ FreshJuice2_pos_54 <dbl> NA, 8141.861, 6760.703, NA, 5957.557, 52588.240, 48…
    $ `40TT_10M2_pos_16` <dbl> NA, 7937.986, 8567.628, NA, 15034.038, 55115.530, 1…
    $ `40TT_10M5_pos_40` <dbl> NA, 5563.556, 8742.438, NA, 17859.370, 69650.164, 1…
    $ `40TT_10M6_pos_72` <dbl> NA, 9136.893, 8224.797, NA, 23304.758, 70638.445, 1…
    $ QC6_pos_43         <dbl> 12426.765, 16833.031, 7340.395, 3263.205, 18441.256…
    $ MEF1_pos_41        <dbl> 46824.566, 40835.254, 8063.132, 28502.540, 19066.36…
    $ QC9_pos_67         <dbl> 13881.810, 14015.707, 9651.587, 4294.461, 20762.291…
    $ MEF_SS6_pos_75     <dbl> 54912.113, 18580.879, 8685.794, 9954.626, 24796.635…
    $ MEF4_pos_29        <dbl> 52074.332, 34421.723, 8003.345, 14098.615, 24177.24…
    $ MEF_SS4_pos_68     <dbl> 61980.793, 22730.098, 7116.988, 7814.761, 23002.707…
    $ MEF_SS3_pos_52     <dbl> 62247.360, 23359.580, 6889.786, 4947.561, 23005.969…
    $ MEF_SS1_pos_56     <dbl> 42746.410, 23413.078, 8209.089, 7750.971, 21103.540…
    $ MEF_SS2_pos_34     <dbl> 61675.145, 20732.110, 6245.571, 10057.129, 21043.80…
    $ MEF3_pos_44        <dbl> 48735.613, 33073.800, 7842.392, 12948.969, 16841.21…
    $ MEF5_pos_50        <dbl> 46009.055, 32630.140, 7295.025, 10213.580, 22322.04…
    $ QC8_pos_59         <dbl> 10166.569, 11754.033, 7449.386, 3190.422, 16414.750…
    $ MEF2_pos_60        <dbl> 48321.387, 37912.586, 7204.479, 15970.995, 23917.12…
    $ QC10_pos_74        <dbl> 12866.921, 14368.647, 8746.157, 3739.938, 21310.840…
    $ QC7_pos_51         <dbl> 14680.529, 16170.804, 11082.852, 1881.733, 15499.25…
    $ QC2_pos_13         <dbl> 12250.746, 14694.995, 5962.542, 3262.584, 15663.376…
    $ QC3_pos_20         <dbl> 13049.109, 17006.514, 6547.161, 3115.725, 17671.970…
    $ QC11_pos_78        <dbl> 14373.564, 14448.783, 9465.725, 4620.273, 19552.660…
    $ QC5_pos_35         <dbl> 16023.473, 12852.819, 6677.819, 3708.668, 16252.766…
    $ QC1_pos_12         <dbl> 14941.828, 11572.707, 7141.102, 4035.967, 15206.507…
    $ QC4_pos_28         <dbl> 12876.858, 13565.415, 8116.300, 3055.272, 15434.733…
    $ QC12_pos_79        <dbl> 14667.275, 12310.267, 7610.508, 3196.210, 23105.748…

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

              sample_name 100.0394_0.4966 100.0748_0.3036 100.112_1.8034
    2  FreshJuice5_pos_73              NA       11359.062       9083.935
    3  FreshJuice3_pos_64              NA        6991.805      10544.096
    4    40TT_10M4_pos_26              NA        9396.625       7078.440
    5  FreshJuice4_pos_38              NA        7958.649       9853.188
    6  FreshJuice6_pos_45              NA        9622.346       9298.826
    7  FreshJuice1_pos_24              NA        7768.952       8173.136
    8    40TT_10M3_pos_46              NA        7989.313       7782.835
    9  FreshJuice2_pos_54              NA        8141.861       6760.703
    10   40TT_10M2_pos_16              NA        7937.986       8567.628
    11   40TT_10M5_pos_40              NA        5563.556       8742.438
       100.1121_1.5856 101.0232_1.6779 101.0233_0.4368 101.0233_1.5452
    2               NA        5578.292        53138.73        5868.525
    3               NA        7092.719        61444.94        5340.145
    4               NA       14991.171        64184.21       10671.049
    5         1311.417        5236.963        49287.32        4020.112
    6               NA        4703.250        47034.38        4744.825
    7               NA        3129.431        51654.85        4754.987
    8               NA       19899.783        67361.87       10135.326
    9               NA        5957.557        52588.24        4812.681
    10              NA       15034.038        55115.53       11640.987
    11              NA       17859.370        69650.16       11810.327
       101.0234_0.8025 102.0372_0.693
    2         16965.26       8268.922
    3         14488.80       6488.537
    4         15152.90       6782.827
    5         17414.89       6841.243
    6         15508.94       8841.724
    7         16636.98       6065.969
    8         15092.41       7913.513
    9         17716.71       6458.659
    10        17292.70       7718.880
    11        18275.40       9692.725

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

    Rows: 108,537
    Columns: 9
    $ sample_name       <chr> "FreshJuice5_pos_73", "FreshJuice5_pos_73", "FreshJu…
    $ short_sample_name <chr> "FJ-5", "FJ-5", "FJ-5", "FJ-5", "FJ-5", "FJ-5", "FJ-…
    $ treatment         <chr> "FJ", "FJ", "FJ", "FJ", "FJ", "FJ", "FJ", "FJ", "FJ"…
    $ run_order         <dbl> 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, …
    $ soluble_solids    <dbl> 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3, 10.3…
    $ mass              <dbl> 0.4268, 0.4268, 0.4268, 0.4268, 0.4268, 0.4268, 0.42…
    $ sample_or_qc      <chr> "Sample", "Sample", "Sample", "Sample", "Sample", "S…
    $ mz_rt             <chr> "100.0394_0.4966", "100.0748_0.3036", "100.112_1.803…
    $ rel_abund         <dbl> NA, 11359.062, 9083.935, NA, 5578.292, 53138.734, 58…

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
      sample_name        short_sample_name treatment run_order soluble_solids  mass
      <chr>              <chr>             <chr>         <dbl>          <dbl> <dbl>
    1 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
    2 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
    3 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
    4 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
    5 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
    6 FreshJuice5_pos_73 FJ-5              FJ               73           10.3 0.427
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
    1 FJ-5        FJ               73 Sample        100. 0.497 100.0394_0…       NA 
    2 FJ-5        FJ               73 Sample        100. 0.304 100.0748_0…   103357.
    3 FJ-5        FJ               73 Sample        100. 1.80  100.112_1.…    82656.
    4 FJ-5        FJ               73 Sample        100. 1.59  100.1121_1…       NA 
    5 FJ-5        FJ               73 Sample        101. 1.68  101.0232_1…    50757.
    6 FJ-5        FJ               73 Sample        101. 0.437 101.0233_0…   483515.

## Data summaries

What mass range do I have?

``` r
range(metab_corrected_long$mz)
```

    [1]  100.0394 1443.3387

What retention time range do I have?

``` r
range(metab_corrected_long$rt)
```

    [1] 0.3004 7.9729

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
       subtitle = "C18 reversed phase, positive ionization mode")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-13-1.png)

Distribution of masses

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>%
  ggplot(aes(x = mz)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  labs(x = "m/z",
       y = "Number of features",
       title = "Distribution of features by mass")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-14-1.png)

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

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-15-1.png)

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
       y = "How many features have that \nmany missing values")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-16-1.png)

How many features have no missing values?

``` r
na_by_feature %>%
  count(missing_values == 0)
```

      missing_values == 0    n
    1               FALSE  481
    2                TRUE 2808

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
       y = "How many samples have that \nmany missing values")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-18-1.png)

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
    1 147.0352_4.5939 TRUE                  21
    2 158.0271_4.5932 TRUE                  21
    3 214.0896_4.5925 TRUE                  21
    4 214.0999_4.5954 TRUE                  21
    5 220.1461_1.7441 TRUE                  21
    6 289.1643_5.498  TRUE                  21

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
    1 FJ-3        TRUE                 223
    2 FJ-5        TRUE                 222
    3 FJ-1        TRUE                 220
    4 FJ-4        TRUE                 220
    5 FJ-6        TRUE                 218
    6 FJ-2        TRUE                 217

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
    1 FJ-5        FJ               73 Sample       100.0394_0.4966    48308.
    2 FJ-5        FJ               73 Sample       100.0748_0.3036   103357.
    3 FJ-5        FJ               73 Sample       100.112_1.8034     82656.
    4 FJ-5        FJ               73 Sample       100.1121_1.5856     5830.
    5 FJ-5        FJ               73 Sample       101.0232_1.6779    50757.
    6 FJ-5        FJ               73 Sample       101.0233_0.4368   483515.

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
       title = "m/z by retention time plot before notame, 2897 features",
       subtitle = "C18 reverse phase, positive ionization mode"))
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-28-1.png)

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

    # A tibble: 33 × 3,290
       sample_name `100.0394_0.4966` `100.0748_0.3036` `100.112_1.8034`
       <chr>                   <dbl>             <dbl>            <dbl>
     1 FJ-5                   48308.           103357.           82656.
     2 FJ-3                   48308.            61879.           93318.
     3 MT-4                   48308.            84937.           63983.
     4 FJ-4                   48308.            70759.           87603.
     5 FJ-6                   48308.            87988.           85029.
     6 FJ-1                   48308.            69969.           73609.
     7 MT-3                   48308.            72474.           70601.
     8 FJ-2                   48308.            73910.           61373.
     9 MT-2                   48308.            71959.           77667.
    10 MT-5                   48308.            50754.           79754.
    # ℹ 23 more rows
    # ℹ 3,286 more variables: `100.1121_1.5856` <dbl>, `101.0232_1.6779` <dbl>,
    #   `101.0233_0.4368` <dbl>, `101.0233_1.5452` <dbl>, `101.0234_0.8025` <dbl>,
    #   `102.0372_0.693` <dbl>, `102.0549_0.5532` <dbl>, `102.055_3.1638` <dbl>,
    #   `102.0914_1.9157` <dbl>, `102.1276_2.1675` <dbl>, `103.039_1.1011` <dbl>,
    #   `103.0391_0.8388` <dbl>, `103.0542_1.6945` <dbl>, `104.0528_0.5517` <dbl>,
    #   `104.0698_0.3039` <dbl>, `104.0706_0.4179` <dbl>, `104.107_5.7303` <dbl>, …

Create df with features.

``` r
features <- metab_corrected_long_imputed %>%
  dplyr::select(mz_rt, mz, rt) %>%
  mutate(across(c(mz, rt), as.numeric)) %>%
  as.data.frame() %>%
  distinct()

glimpse(features)
```

    Rows: 3,289
    Columns: 3
    $ mz_rt <chr> "100.0394_0.4966", "100.0748_0.3036", "100.112_1.8034", "100.112…
    $ mz    <dbl> 100.0394, 100.0748, 100.1120, 100.1121, 101.0232, 101.0233, 101.…
    $ rt    <dbl> 0.4966, 0.3036, 1.8034, 1.5856, 1.6779, 0.4368, 1.5452, 0.8025, …

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
    [1] 3000
    [1] 3100
    [1] 3200

``` r
head(connection)
```

                    x               y       cor rt_diff mz_diff
    1 100.0394_0.4966 133.0497_0.4822 0.9678008 -0.0144 33.0103
    2 100.0394_0.4966 160.1332_0.4997 0.9715157  0.0031 60.0938
    3 100.0748_0.3036  159.075_0.3034 0.9633314 -0.0002 59.0002
    4 100.0748_0.3036 161.0911_0.3188 0.9549384  0.0152 61.0163
    5 100.0748_0.3036 176.1012_0.3036 0.9656791  0.0000 76.0264
    6 101.0232_1.6779 143.0339_1.6772 0.9569108 -0.0007 42.0107

### Clustering

Now that we have found all of the features that are connected based on
the parameters we have set, we now need to find clusters.

``` r
clusters <- find_clusters(connections = connection, 
                          d_thresh = 0.8)
```

    234 components found

    Component 100 / 234 
    Component 200 / 234 
    170 components found

    Component 100 / 170 
    41 components found

    18 components found

    7 components found

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
          "features_clustered_pos_MEF.csv")
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

    [1] 3289

``` r
# how many features got removed during clustering?
nrow(metab) - nrow(cluster_features)
```

    [1] 822

``` r
# what percentage of the original features were removed?
((nrow(metab) - nrow(cluster_features))/nrow(metab)) * 100
```

    [1] 24.9924

Reduce our dataset to include only our new clusters. `cluster_data`
contains only the retained clusters, while `cluster_features` tells you
also which features are a part of each cluster.

``` r
# combined metadata_plus with cluster_features

metab_imputed_clustered_wide <- left_join(metab_corrected_wide_imputed[,1:4], cluster_data,
                                          by = "sample_name") 

dim(metab_imputed_clustered_wide) # we have 2474 features since 4 metadata columns
```

    [1]   33 2471

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

Write out the final clustered dataset

``` r
write_csv(metab_imputed_clustered_long,
  "features_clustered_long__pos_MEF.csv")
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
       title = "m/z by retention time plot after notame, 2219 features",
       subtitle = "C18 reverse phase, positive ionization mode"))
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-39-1.png)

``` r
before_notame / after_notame
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-40-1.png)

## Assessing data quality

Let’s make sure that our data is of good quality.

### Untransformed data

First we are going to convert the type of some of the columns to match
what we want (e.g., run order converted to numeric, treatment to factor)

``` r
tibble(metab_imputed_clustered_long)
```

    # A tibble: 81,411 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <chr>         <dbl> <chr>        <dbl> <dbl> <chr>          <dbl>
     1 FJ-5        FJ               73 Sample        100. 1.80  100.112_1…    82656.
     2 FJ-5        FJ               73 Sample        100. 1.59  100.1121_…     5830.
     3 FJ-5        FJ               73 Sample        101. 0.437 101.0233_…   483515.
     4 FJ-5        FJ               73 Sample        101. 1.55  101.0233_…    53398.
     5 FJ-5        FJ               73 Sample        101. 0.802 101.0234_…   154369.
     6 FJ-5        FJ               73 Sample        102. 0.693 102.0372_…    75240.
     7 FJ-5        FJ               73 Sample        102. 1.92  102.0914_…    20643.
     8 FJ-5        FJ               73 Sample        102. 2.17  102.1276_…   697850.
     9 FJ-5        FJ               73 Sample        103. 1.10  103.039_1…   203268.
    10 FJ-5        FJ               73 Sample        103. 0.839 103.0391_…   265538.
    # ℹ 81,401 more rows

``` r
# make run_order numeric
metab_imputed_clustered_long$run_order <- as.numeric(metab_imputed_clustered_long$run_order)

# make treatment and sample_or_qc a factor (i.e., categorical)
metab_imputed_clustered_long$treatment <- as.factor(metab_imputed_clustered_long$treatment)
metab_imputed_clustered_long$sample_or_qc <- as.factor(metab_imputed_clustered_long$sample_or_qc)

# did it work?
tibble(metab_imputed_clustered_long)
```

    # A tibble: 81,411 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <fct>         <dbl> <fct>        <dbl> <dbl> <chr>          <dbl>
     1 FJ-5        FJ               73 Sample        100. 1.80  100.112_1…    82656.
     2 FJ-5        FJ               73 Sample        100. 1.59  100.1121_…     5830.
     3 FJ-5        FJ               73 Sample        101. 0.437 101.0233_…   483515.
     4 FJ-5        FJ               73 Sample        101. 1.55  101.0233_…    53398.
     5 FJ-5        FJ               73 Sample        101. 0.802 101.0234_…   154369.
     6 FJ-5        FJ               73 Sample        102. 0.693 102.0372_…    75240.
     7 FJ-5        FJ               73 Sample        102. 1.92  102.0914_…    20643.
     8 FJ-5        FJ               73 Sample        102. 2.17  102.1276_…   697850.
     9 FJ-5        FJ               73 Sample        103. 1.10  103.039_1…   203268.
    10 FJ-5        FJ               73 Sample        103. 0.839 103.0391_…   265538.
    # ℹ 81,401 more rows

Let’s make a boxplot to see how the metabolite abundance looks across
different samples.

``` r
metab_imputed_clustered_long %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is unscaled",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-42-1.png)

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
metab_imputed_clustered_long_log2 %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is log2 transformed",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-44-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
metab_imputed_clustered_long_log2 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is log2 transformed",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-45-1.png)

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
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-47-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
metab_imputed_clustered_long_log10 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-48-1.png)

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
       sample_name treatment run_order sample_or_qc `100.112_1.8034`[,1]
       <chr>       <fct>         <dbl> <fct>                       <dbl>
     1 FJ-5        FJ               73 Sample                      0.547
     2 FJ-3        FJ               64 Sample                      1.50 
     3 MT-4        MT               26 Sample                     -1.12 
     4 FJ-4        FJ               38 Sample                      0.989
     5 FJ-6        FJ               45 Sample                      0.759
     6 FJ-1        FJ               24 Sample                     -0.262
     7 MT-3        MT               46 Sample                     -0.530
     8 FJ-2        FJ               54 Sample                     -1.35 
     9 MT-2        MT               16 Sample                      0.101
    10 MT-5        MT               40 Sample                      0.287
    # ℹ 5 more variables: `100.1121_1.5856` <dbl[,1]>, `101.0233_0.4368` <dbl[,1]>,
    #   `101.0233_1.5452` <dbl[,1]>, `101.0234_0.8025` <dbl[,1]>,
    #   `102.0372_0.693` <dbl[,1]>

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
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is autoscaled",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-50-1.png)

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
  labs(title = "LC-MS (+) Feature Abundances by Sample",
       subtitle = "Data is Pareto scaled",
       y = "Relative abundance")
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-51-1.png)

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
          file = "features_clustered_wide_log2_pos_MEF.csv")
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
                               PC1     PC2     PC3     PC4    PC5     PC6     PC7
    Standard deviation     24.9429 12.5424 7.28981 5.45165 3.0786 2.79099 2.66254
    Proportion of Variance  0.6436  0.1627 0.05497 0.03075 0.0098 0.00806 0.00733
    Cumulative Proportion   0.6436  0.8063 0.86130 0.89205 0.9019 0.90991 0.91724
                               PC8     PC9    PC10    PC11    PC12    PC13    PC14
    Standard deviation     2.49310 2.33875 2.16049 2.04880 2.03239 1.96505 1.91872
    Proportion of Variance 0.00643 0.00566 0.00483 0.00434 0.00427 0.00399 0.00381
    Cumulative Proportion  0.92367 0.92933 0.93416 0.93850 0.94278 0.94677 0.95058
                              PC15    PC16    PC17    PC18    PC19   PC20    PC21
    Standard deviation     1.88115 1.85052 1.79732 1.78167 1.75553 1.7039 1.69102
    Proportion of Variance 0.00366 0.00354 0.00334 0.00328 0.00319 0.0030 0.00296
    Cumulative Proportion  0.95424 0.95778 0.96112 0.96441 0.96760 0.9706 0.97356
                              PC22   PC23    PC24    PC25    PC26    PC27   PC28
    Standard deviation     1.68025 1.6443 1.61948 1.56327 1.52036 1.51159 1.4924
    Proportion of Variance 0.00292 0.0028 0.00271 0.00253 0.00239 0.00236 0.0023
    Cumulative Proportion  0.97648 0.9793 0.98199 0.98452 0.98691 0.98927 0.9916
                              PC29    PC30   PC31    PC32      PC33
    Standard deviation     1.46594 1.44018 1.4254 1.37473 4.833e-14
    Proportion of Variance 0.00222 0.00215 0.0021 0.00196 0.000e+00
    Cumulative Proportion  0.99380 0.99594 0.9980 1.00000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_qc <- summary(pca_qc)$importance %>%
  as.data.frame()

head(importance_qc)
```

                                PC1      PC2     PC3      PC4      PC5      PC6
    Standard deviation     24.94289 12.54241 7.28981 5.451649 3.078599 2.790994
    Proportion of Variance  0.64359  0.16273 0.05497 0.030750 0.009800 0.008060
    Cumulative Proportion   0.64359  0.80633 0.86130 0.892050 0.901850 0.909910
                                PC7      PC8      PC9    PC10     PC11     PC12
    Standard deviation     2.662543 2.493101 2.338747 2.16049 2.048803 2.032389
    Proportion of Variance 0.007330 0.006430 0.005660 0.00483 0.004340 0.004270
    Cumulative Proportion  0.917240 0.923670 0.929330 0.93416 0.938500 0.942780
                               PC13    PC14     PC15    PC16     PC17     PC18
    Standard deviation     1.965052 1.91872 1.881152 1.85052 1.797324 1.781674
    Proportion of Variance 0.003990 0.00381 0.003660 0.00354 0.003340 0.003280
    Cumulative Proportion  0.946770 0.95058 0.954240 0.95778 0.961120 0.964410
                               PC19     PC20     PC21     PC22     PC23     PC24
    Standard deviation     1.755526 1.703891 1.691021 1.680248 1.644279 1.619477
    Proportion of Variance 0.003190 0.003000 0.002960 0.002920 0.002800 0.002710
    Cumulative Proportion  0.967600 0.970600 0.973560 0.976480 0.979270 0.981990
                               PC25     PC26    PC27     PC28     PC29     PC30
    Standard deviation     1.563271 1.520357 1.51159 1.492435 1.465943 1.440182
    Proportion of Variance 0.002530 0.002390 0.00236 0.002300 0.002220 0.002150
    Cumulative Proportion  0.984520 0.986910 0.98927 0.991570 0.993800 0.995940
                               PC31     PC32         PC33
    Standard deviation     1.425379 1.374729 4.833092e-14
    Proportion of Variance 0.002100 0.001960 0.000000e+00
    Cumulative Proportion  0.998040 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-55-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-56-1.png)

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
       title = "PCA Scores Plot Colored by Treatment, positive Mode"))
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-58-1.png)

Then make your scores plot ineractive so you can see who is who.

``` r
ggplotly(scores_qc_plot, tooltip = "text")
```

<div class="plotly html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-4753da52ceccf21a564c" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-4753da52ceccf21a564c">{"x":{"data":[{"x":[-33.6828366822697,47.2346378049496],"y":[0,0],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,0],"y":[-23.5954956469483,13.1951786827926],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-30.0047696601233,-29.4839548772305,-29.5206024622424,-29.8537622720394,-29.284807387557,-29.1209634089234],"y":[-20.5666966316187,-21.1904229106161,-21.9231922683238,-21.3731129998547,-21.8420844999807,-20.6067433815097],"text":["Sample: FJ-5,<br />Treatment: FJ","Sample: FJ-3,<br />Treatment: FJ","Sample: FJ-4,<br />Treatment: FJ","Sample: FJ-6,<br />Treatment: FJ","Sample: FJ-1,<br />Treatment: FJ","Sample: FJ-2,<br />Treatment: FJ"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"FJ","legendgroup":"FJ","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[43.5565707828032,27.7395021593634,28.6894968089495,25.0533585399433,32.2533601103378],"y":[-9.78707545178355,-1.93195721014533,-4.23260062428322,-0.436430142233261,-4.46765332881572],"text":["Sample: MEF-1,<br />Treatment: MEF","Sample: MEF-4,<br />Treatment: MEF","Sample: MEF-3,<br />Treatment: MEF","Sample: MEF-5,<br />Treatment: MEF","Sample: MEF-2,<br />Treatment: MEF"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(163,165,0,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF","legendgroup":"MEF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[39.0190189522171,37.7263843495787,38.0967403687244,37.4549994498309,37.8643583515105],"y":[-7.29783370205252,-10.4393008615139,-6.61969673645315,-3.84155137571015,-4.68071174427266],"text":["Sample: MEF+SS-6,<br />Treatment: MEF+SS","Sample: MEF+SS-4,<br />Treatment: MEF+SS","Sample: MEF+SS-3,<br />Treatment: MEF+SS","Sample: MEF+SS-1,<br />Treatment: MEF+SS","Sample: MEF+SS-2,<br />Treatment: MEF+SS"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,191,125,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF+SS","legendgroup":"MEF+SS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-18.652508199223,-19.1357360550046,-18.9552157547583,-19.1010276550584,-19.2350818243428],"y":[10.3189032140241,10.7850824430966,8.93981591773362,9.43303343520874,9.54550358411694],"text":["Sample: MT-4,<br />Treatment: MT","Sample: MT-3,<br />Treatment: MT","Sample: MT-2,<br />Treatment: MT","Sample: MT-5,<br />Treatment: MT","Sample: MT-6,<br />Treatment: MT"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,176,246,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MT","legendgroup":"MT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-6.32995773897477,-6.27186090192305,-6.33544303306926,-6.6544665549101,-6.23603909558528,-5.94800937185694,-5.91186521893512,-6.24626813693254,-6.02272552284894,-5.81215512012392,-6.43407746990741,-6.90249215168802],"y":[11.0874387013455,10.8407877331436,10.7311804896374,11.3221260195483,11.522875304168,10.771887839288,11.469908411809,11.3559785484292,11.1088831595638,11.2050370667357,11.2074037733951,9.59121822792333],"text":["Sample: QC-6,<br />Treatment: QC","Sample: QC-9,<br />Treatment: QC","Sample: QC-8,<br />Treatment: QC","Sample: QC-10,<br />Treatment: QC","Sample: QC-7,<br />Treatment: QC","Sample: QC-2,<br />Treatment: QC","Sample: QC-3,<br />Treatment: QC","Sample: QC-11,<br />Treatment: QC","Sample: QC-5,<br />Treatment: QC","Sample: QC-1,<br />Treatment: QC","Sample: QC-4,<br />Treatment: QC","Sample: QC-12,<br />Treatment: QC"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(231,107,243,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"QC","legendgroup":"QC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":43.7625570776256,"r":7.30593607305936,"b":40.1826484018265,"l":43.1050228310502},"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"title":{"text":"PCA Scores Plot Colored by Treatment, positive Mode","font":{"color":"rgba(0,0,0,1)","family":"","size":17.5342465753425},"x":0,"xref":"paper"},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-33.6828366822697,47.2346378049496],"tickmode":"array","ticktext":["-20","0","20","40"],"tickvals":[-20,0,20,40],"categoryorder":"array","categoryarray":["-20","0","20","40"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"PC1: 64.4%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-23.5954956469483,13.1951786827926],"tickmode":"array","ticktext":["-20","-10","0","10"],"tickvals":[-20,-10,0,10],"categoryorder":"array","categoryarray":["-20","-10","0","10"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"PC2: 16.3%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"title":{"text":"treatment","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"db344bc33e60":{"yintercept":{},"type":"scatter"},"db3431096065":{"xintercept":{}},"db3425804ae":{"x":{},"y":{},"fill":{},"text":{}}},"cur_data":"db344bc33e60","visdat":{"db344bc33e60":["function (y) ","x"],"db3431096065":["function (y) ","x"],"db3425804ae":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-60-1.png)

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
                               PC1     PC2     PC3     PC4     PC5     PC6    PC7
    Standard deviation     31.0601 13.0388 7.25584 3.95986 3.37767 3.33825 3.0261
    Proportion of Variance  0.7377  0.1300 0.04026 0.01199 0.00872 0.00852 0.0070
    Cumulative Proportion   0.7377  0.8677 0.90798 0.91997 0.92869 0.93722 0.9442
                               PC8     PC9    PC10    PC11    PC12    PC13    PC14
    Standard deviation     2.85719 2.67523 2.57421 2.53639 2.43134 2.39312 2.32847
    Proportion of Variance 0.00624 0.00547 0.00507 0.00492 0.00452 0.00438 0.00415
    Cumulative Proportion  0.95046 0.95593 0.96100 0.96592 0.97044 0.97482 0.97897
                            PC15    PC16    PC17    PC18    PC19   PC20      PC21
    Standard deviation     2.288 2.25106 2.18642 2.12734 2.05730 1.9145 4.727e-14
    Proportion of Variance 0.004 0.00387 0.00366 0.00346 0.00324 0.0028 0.000e+00
    Cumulative Proportion  0.983 0.98684 0.99050 0.99396 0.99720 1.0000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_noqc <- summary(pca_noqc)$importance %>%
  as.data.frame()

head(importance_noqc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     31.06014 13.03878 7.255844 3.959861 3.377674 3.338248
    Proportion of Variance  0.73772  0.13000 0.040260 0.011990 0.008720 0.008520
    Cumulative Proportion   0.73772  0.86772 0.907980 0.919970 0.928690 0.937220
                                PC7     PC8      PC9    PC10     PC11     PC12
    Standard deviation     3.026114 2.85719 2.675232 2.57421 2.536388 2.431342
    Proportion of Variance 0.007000 0.00624 0.005470 0.00507 0.004920 0.004520
    Cumulative Proportion  0.944220 0.95046 0.955930 0.96100 0.965920 0.970440
                               PC13     PC14     PC15     PC16     PC17     PC18
    Standard deviation     2.393119 2.328466 2.288089 2.251056 2.186418 2.127341
    Proportion of Variance 0.004380 0.004150 0.004000 0.003870 0.003660 0.003460
    Cumulative Proportion  0.974820 0.978970 0.982970 0.986840 0.990500 0.993960
                               PC19     PC20         PC21
    Standard deviation     2.057302 1.914535 4.726953e-14
    Proportion of Variance 0.003240 0.002800 0.000000e+00
    Cumulative Proportion  0.997200 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-63-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-64-1.png)

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
       title = "PCA Scores Plot Colored by Treatment, Positive Mode"))
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-66-1.png)

Then make your scores plot ineractive so you can see who is who.

``` r
ggplotly(scores_noqc_plot, tooltip = "text")
```

<div class="plotly html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-1319471b434b3d1b5eba" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-1319471b434b3d1b5eba">{"x":{"data":[{"x":[-38.0399378870177,43.3644171899089],"y":[0,0],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,0],"y":[-16.8985483781534,23.333653438476],"text":"","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-34.3397399289755,-33.8839372316696,-33.9713886681734,-34.2474021709413,-33.6903918699513,-33.4793289622314],"y":[-13.5126845759174,-13.8743142679104,-14.5623629689983,-14.3149874525551,-15.0698119319429,-13.3747317087594],"text":["Sample: FJ-5,<br />Treatment: FJ","Sample: FJ-3,<br />Treatment: FJ","Sample: FJ-4,<br />Treatment: FJ","Sample: FJ-6,<br />Treatment: FJ","Sample: FJ-1,<br />Treatment: FJ","Sample: FJ-2,<br />Treatment: FJ"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"FJ","legendgroup":"FJ","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[39.6642192318667,24.5277563802354,25.3181186835047,21.9475268701639,28.8340021871359],"y":[-6.25629292173874,0.474572191848113,-1.7248418034217,2.1416535111945,-1.73991839708234],"text":["Sample: MEF-1,<br />Treatment: MEF","Sample: MEF-4,<br />Treatment: MEF","Sample: MEF-3,<br />Treatment: MEF","Sample: MEF-5,<br />Treatment: MEF","Sample: MEF-2,<br />Treatment: MEF"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(124,174,0,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF","legendgroup":"MEF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[35.2943261239168,33.8001815437814,34.3698970253598,33.9938994289986,34.3181929370737],"y":[-2.87801309912679,-5.97681915643207,-1.56363035076108,0.0932114097088413,-0.391832140694809],"text":["Sample: MEF+SS-6,<br />Treatment: MEF+SS","Sample: MEF+SS-4,<br />Treatment: MEF+SS","Sample: MEF+SS-3,<br />Treatment: MEF+SS","Sample: MEF+SS-1,<br />Treatment: MEF+SS","Sample: MEF+SS-2,<br />Treatment: MEF+SS"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MEF+SS","legendgroup":"MEF+SS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-21.3315056349958,-21.7280234571622,-21.6665522411037,-21.8000656368296,-21.9297846100034],"y":[21.5049169922656,21.2532539085348,19.3609825732005,20.1234569687455,20.2881932198432],"text":["Sample: MT-4,<br />Treatment: MT","Sample: MT-3,<br />Treatment: MT","Sample: MT-2,<br />Treatment: MT","Sample: MT-5,<br />Treatment: MT","Sample: MT-6,<br />Treatment: MT"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(199,124,255,1)","opacity":1,"size":5.66929133858268,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)"}},"hoveron":"points","name":"MT","legendgroup":"MT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":43.7625570776256,"r":7.30593607305936,"b":40.1826484018265,"l":43.1050228310502},"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"title":{"text":"PCA Scores Plot Colored by Treatment, Positive Mode","font":{"color":"rgba(0,0,0,1)","family":"","size":17.5342465753425},"x":0,"xref":"paper"},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-38.0399378870177,43.3644171899089],"tickmode":"array","ticktext":["-20","0","20","40"],"tickvals":[-20,0,20,40],"categoryorder":"array","categoryarray":["-20","0","20","40"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"PC1: 73.8%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-16.8985483781534,23.333653438476],"tickmode":"array","ticktext":["-10","0","10","20"],"tickvals":[-10,0,10,20],"categoryorder":"array","categoryarray":["-10","0","10","20"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"PC2: 13%","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"title":{"text":"treatment","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"db3476874b4c":{"yintercept":{},"type":"scatter"},"db3449ba285d":{"xintercept":{}},"db342b7ecb5":{"x":{},"y":{},"fill":{},"text":{}}},"cur_data":"db3476874b4c","visdat":{"db3476874b4c":["function (y) ","x"],"db3449ba285d":["function (y) ","x"],"db342b7ecb5":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_noqc)
```

![](Metabolomics---Pos-Data---MEF---Data-Quality_files/figure-commonmark/unnamed-chunk-68-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.
