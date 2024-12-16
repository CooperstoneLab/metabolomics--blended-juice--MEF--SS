# Filtering - Untargeted Metabolomics Data - Positive Mode
Giovana Domeneghini Mercali
2023-08-30

# Loading libraries

``` r
library(tidyverse)
library(janitor) # if you want to clean_names()
```

Once you get deconvoluted data from MZmine or similar programs, you need
to wrangle your data in such a way that you can conduct your analysis on
it.

# Read in data

First we want to read in our raw data. The code here is to read in data
directly from MZmine, so if you are using a different program for
deconvolution, you might need to make some light adjustments.

``` r
metadata_juice <- read_csv(file = "Feature_list_all_data_pos.csv",
                      col_names = TRUE, # has headers
                      na = "0")
```

    Rows: 6546 Columns: 82
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    dbl (82): row ID, row m/z, row retention time, ConditioningQC_pos_3.mzML Pea...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# How many features do we have?
dim(metadata_juice)
```

    [1] 6546   82

``` r
# Look at beginning of the dataset

metadata_juice[1:8, 1:8]
```

    # A tibble: 8 × 8
      `row ID` `row m/z` `row retention time` ConditioningQC_pos_3.mzML Peak heigh…¹
         <dbl>     <dbl>                <dbl>                                  <dbl>
    1      906      100.                0.497                                 13687.
    2      246      100.                0.304                                 10000.
    3     2700      100.                1.74                                  14797.
    4     2770      100.                1.80                                   7451.
    5     2472      100.                1.59                                     NA 
    6     2842      101.                1.81                                   4433.
    7     2323      101.                1.44                                  11262.
    8     2109      101.                1.08                                  15022.
    # ℹ abbreviated name: ¹​`ConditioningQC_pos_3.mzML Peak height`
    # ℹ 4 more variables: `ConditioningQC_pos_10.mzML Peak height` <dbl>,
    #   `ConditioningQC_pos_5.mzML Peak height` <dbl>,
    #   `40TT2_pos_22.mzML Peak height` <dbl>,
    #   `ConditioningQC_pos_9.mzML Peak height` <dbl>

``` r
# str(metadata_juice)
# names(metadata_juice)
#length(metadata_juice)
```

Note there is no metadata included in this file. Just m/z, retention
time, and a column for each sample, where values are peak heights. We
are using peak height instead of peak area because it is less dependent
on bad peak shape which you get sometimes with Metabolomics.

# RT filter

You might have deconvoluted more data than you plan to use in your
analysis. For example, you may want to exclude the first bit and last
bit of your run, since you do not expect to have good reproducibility in
those areas.

Here, we are filtering to only include features that elute between
0.3-8.0 min of this 10 min run.

``` r
metadata_juice_RTfilt <- metadata_juice %>%
  filter(between(`row retention time`, 0.3, 8.0))

# Did it work?
range(metadata_juice_RTfilt$`row retention time`)
```

    [1] 0.3004379 7.9969680

``` r
# How many features do we have now?
dim(metadata_juice_RTfilt)
```

    [1] 5783   82

# Cleaning up data

## Create mz_rt

This creates a unique identifier for each feature using its
mass-to-charge ratio (m/z) and retention time (RT).

``` r
MZ_RT <- metadata_juice_RTfilt %>%
  mutate(mz = round(metadata_juice_RTfilt$`row m/z`, digits = 4), # decrease number of decimals for m/z & rt
         rt = round(metadata_juice_RTfilt$`row retention time`, digits = 4),
         .before = 1,
         .keep = "unused") %>%
  unite(mz_rt, c(mz, rt), remove = TRUE) %>% # Combine m/z & rt with _ in between
  rename(mz = `row m/z`,
         rt = `row retention time`,
         row_ID = `row ID`) %>% # rename to be more R friendly
  select(row_ID, everything()) # reorder and move row_ID to front

MZ_RT[1:8, 1:8]
```

    # A tibble: 8 × 8
      row_ID mz_rt            mz    rt ConditioningQC_pos_3…¹ ConditioningQC_pos_1…²
       <dbl> <chr>         <dbl> <dbl>                  <dbl>                  <dbl>
    1    906 100.0394_0.4…  100. 0.497                 13687.                 12728.
    2    246 100.0748_0.3…  100. 0.304                 10000.                 13485.
    3   2700 100.0757_1.7…  100. 1.74                  14797.                 13501.
    4   2770 100.112_1.80…  100. 1.80                   7451.                  8237.
    5   2472 100.1121_1.5…  100. 1.59                     NA                   2431.
    6   2842 101.0072_1.8…  101. 1.81                   4433.                  6216.
    7   2323 101.0072_1.4…  101. 1.44                  11262.                 13175.
    8   2109 101.0072_1.0…  101. 1.08                  15022.                 15879.
    # ℹ abbreviated names: ¹​`ConditioningQC_pos_3.mzML Peak height`,
    #   ²​`ConditioningQC_pos_10.mzML Peak height`
    # ℹ 2 more variables: `ConditioningQC_pos_5.mzML Peak height` <dbl>,
    #   `40TT2_pos_22.mzML Peak height` <dbl>

## Clean up file names

We are using `gsub()` to replace strings (i.e. characters) in our sample
names. Here is some useful info about
[gsub](http://www.endmemo.com/r/gsub.php), and two different tutorials
[here](https://www.youtube.com/watch?v=4ZokHoF99DY) and
[here](https://www.youtube.com/watch?v=r4Sh7H6wzPA). You will likely
need to change this code to suit your purposes.

``` r
# remove stuff from the end of file names, ".mzML Peak height"
newcolumnnames <- gsub(".mzML.*","", colnames(MZ_RT))
colnames(MZ_RT) <- newcolumnnames
```

What are our sample names?

``` r
colnames(MZ_RT)
```

     [1] "row_ID"                "mz_rt"                 "mz"                   
     [4] "rt"                    "ConditioningQC_pos_3"  "ConditioningQC_pos_10"
     [7] "ConditioningQC_pos_5"  "40TT2_pos_22"          "ConditioningQC_pos_9" 
    [10] "FreshJuice5_pos_73"    "40TT4_pos_49"          "FreshJuice3_pos_64"   
    [13] "70TT3_pos_14"          "FreshJuice0_pos_17"    "40TT6_pos_69"         
    [16] "40TT_10M4_pos_26"      "FreshJuice4_pos_38"    "40TT1_pos_53"         
    [19] "ConditioningQC_pos_4"  "ConditioningQC_pos_11" "40TT5_pos_61"         
    [22] "40TT3_pos_39"          "70TT2_pos_25"          "40TT_10M1_pos_18"     
    [25] "FreshJuice6_pos_45"    "FreshJuice1_pos_24"    "40TT_10M3_pos_46"     
    [28] "70TT7_pos_76"          "FreshJuice2_pos_54"    "70TT4_pos_55"         
    [31] "40TT_10M2_pos_16"      "70TT5_pos_30"          "70TT6_pos_32"         
    [34] "40TT_10M5_pos_40"      "40TT_10M6_pos_72"      "70TT1_pos_70"         
    [37] "ProcessBlank_pos_8"    "SolventBlank3_pos_80"  "SolventBlank4_pos_83" 
    [40] "ProcessBlank2_pos_77"  "SolventBlank2_pos_2"   "SolventBlank1_pos_1"  
    [43] "QC6_pos_43"            "MEF1_pos_41"           "QC9_pos_67"           
    [46] "MEF_SS6_pos_75"        "MEF_SS5_pos_36"        "MEF4_pos_29"          
    [49] "MEF_SS4_pos_68"        "MEF_SS3_pos_52"        "MEF_SS1_pos_56"       
    [52] "MEF_SS2_pos_34"        "HPP1_pos_33"           "HPP4_pos_21"          
    [55] "HPP6_pos_31"           "MEF3_pos_44"           "MEF5_pos_50"          
    [58] "QC8_pos_59"            "MEF2_pos_60"           "HPP3_pos_27"          
    [61] "HPP2_pos_15"           "QC10_pos_74"           "QC7_pos_51"           
    [64] "HPP5_pos_71"           "QC2_pos_13"            "QC3_pos_20"           
    [67] "QC11_pos_78"           "QC5_pos_35"            "US1_pos_23"           
    [70] "UST2_pos_19"           "QC1_pos_12"            "UST1_pos_37"          
    [73] "US3_pos_47"            "QC4_pos_28"            "UST5_pos_63"          
    [76] "US5_pos_42"            "QC12_pos_79"           "US2_pos_48"           
    [79] "UST6_pos_62"           "US4_pos_65"            "UST3_pos_57"          
    [82] "US6_pos_66"            "UST4_pos_58"          

## Remove unwanted samples

Removing conditioning QCs and solvent blanks.

``` r
MZ_RT <- MZ_RT %>%
  select(-contains("Conditioning"), -contains("Solvent"))

colnames(MZ_RT)
```

     [1] "row_ID"               "mz_rt"                "mz"                  
     [4] "rt"                   "40TT2_pos_22"         "FreshJuice5_pos_73"  
     [7] "40TT4_pos_49"         "FreshJuice3_pos_64"   "70TT3_pos_14"        
    [10] "FreshJuice0_pos_17"   "40TT6_pos_69"         "40TT_10M4_pos_26"    
    [13] "FreshJuice4_pos_38"   "40TT1_pos_53"         "40TT5_pos_61"        
    [16] "40TT3_pos_39"         "70TT2_pos_25"         "40TT_10M1_pos_18"    
    [19] "FreshJuice6_pos_45"   "FreshJuice1_pos_24"   "40TT_10M3_pos_46"    
    [22] "70TT7_pos_76"         "FreshJuice2_pos_54"   "70TT4_pos_55"        
    [25] "40TT_10M2_pos_16"     "70TT5_pos_30"         "70TT6_pos_32"        
    [28] "40TT_10M5_pos_40"     "40TT_10M6_pos_72"     "70TT1_pos_70"        
    [31] "ProcessBlank_pos_8"   "ProcessBlank2_pos_77" "QC6_pos_43"          
    [34] "MEF1_pos_41"          "QC9_pos_67"           "MEF_SS6_pos_75"      
    [37] "MEF_SS5_pos_36"       "MEF4_pos_29"          "MEF_SS4_pos_68"      
    [40] "MEF_SS3_pos_52"       "MEF_SS1_pos_56"       "MEF_SS2_pos_34"      
    [43] "HPP1_pos_33"          "HPP4_pos_21"          "HPP6_pos_31"         
    [46] "MEF3_pos_44"          "MEF5_pos_50"          "QC8_pos_59"          
    [49] "MEF2_pos_60"          "HPP3_pos_27"          "HPP2_pos_15"         
    [52] "QC10_pos_74"          "QC7_pos_51"           "HPP5_pos_71"         
    [55] "QC2_pos_13"           "QC3_pos_20"           "QC11_pos_78"         
    [58] "QC5_pos_35"           "US1_pos_23"           "UST2_pos_19"         
    [61] "QC1_pos_12"           "UST1_pos_37"          "US3_pos_47"          
    [64] "QC4_pos_28"           "UST5_pos_63"          "US5_pos_42"          
    [67] "QC12_pos_79"          "US2_pos_48"           "UST6_pos_62"         
    [70] "US4_pos_65"           "UST3_pos_57"          "US6_pos_66"          
    [73] "UST4_pos_58"         

# Start filtering

### CV function

Since base R does not have a function to calculate coefficient of
variance, let’s write one.

``` r
cv <- function(x){
        (sd(x)/mean(x))
}
```

## Counting QCs

Subset QCs and filter features to keep only those that are present in
100% of QCs. You could change this parameter based on your data.

``` r
# check dimensions of current df
dim(MZ_RT)
```

    [1] 5783   73

``` r
MZ_RT_QCs <- MZ_RT %>%
  select(mz_rt, contains("QC")) %>% # select QCs
  filter(rowSums(is.na(.)) <= 1) # remove rows that have 1 or more NAs
```

``` r
# check dimensions of QCs filtered df
dim(MZ_RT_QCs)
```

    [1] 5632   13

``` r
# how many features got removed with this filtering?
nrow(MZ_RT) - nrow(MZ_RT_QCs)
```

    [1] 151

## Filter on QC CV

Here we are removing features that have a CV of more than 30% in the
QCs. The rationale is that if a feature cannot be reproducibly measured
in samples that are all the same, it should not be included in our
analysis.

``` r
# calculate CV row-wise (1 means row-wise)
QC_CV <- apply(MZ_RT_QCs[, 2:ncol(MZ_RT_QCs)], 1, cv)


# bind the CV vector back to the QC df
MZ_RT_QCs_CV <- cbind(MZ_RT_QCs, QC_CV)

# filter for keeping features with QC_CV <= 0.30 (or 30%)
MZ_RT_QCs_CVfilt <- MZ_RT_QCs_CV %>%
  filter(QC_CV <= 0.30)
```

How many features did I remove with this CV filtering?

``` r
nrow(MZ_RT_QCs) - nrow(MZ_RT_QCs_CVfilt)
```

    [1] 226

## Merge back the rest of the data

MZ_RT_QCs_CVfilt only contains the QCs, We want to keep only the rows
that are present in this df, and then merge back all of the other
samples present in MZ_RT. We will do this by creating a vector that has
the mz_rt features we want to keep, and then using `filter()` and `%in%`
to keep only features that are a part of this list.

``` r
dim(MZ_RT_QCs_CVfilt)
```

    [1] 5406   14

``` r
dim(MZ_RT)
```

    [1] 5783   73

``` r
# make a character vector of the mz_rt features we want to keep
# i.e., the ones that passed our previous filtering steps

features_to_keep <- as.character(MZ_RT_QCs_CVfilt$mz_rt)

MZ_RT_filt <- MZ_RT %>%
  filter(mz_rt %in% features_to_keep)

dim(MZ_RT_filt)
```

    [1] 5406   73

You should have the same number of features in MZ_RT_QCs_CVfilt as you
do in your new filtered df MZ_RT_filt.

``` r
all.equal(nrow(MZ_RT_QCs_CVfilt), nrow(MZ_RT_filt))
```

    [1] TRUE

## Process blanks

We want to remove features that are present in our process blanks as
they are not coming from compounds present in our samples. In this
dataset, the sample (there is only two, typically you would have at
least 3 process blanks to average) representing this process blank (a
sample that includes all the extraction materials, minus the sample,
here the juice was replaced by mass with water) has “PB” in the sample
name.

``` r
# grab the name of the column/sample that is the process blank
grep("ProcessBlank", colnames(MZ_RT_filt), value = TRUE)
```

    [1] "ProcessBlank_pos_8"   "ProcessBlank2_pos_77"

Calculate the average value across the QCs, then remove features that
are not at least 10x higher in the QCs than in the process blank. To do
this we will use
[`apply()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/apply).

`apply(X, MARGIN, FUN,...)` where X is your df, MARGIN is 1 for
row-wise, and 2 for col-wise, and FUN is your function

``` r
# pull avg peak height across QCs
avg_height_QC <- apply(MZ_RT_QCs_CVfilt[, 2:ncol(MZ_RT_QCs_CVfilt)], 1, mean)

# bind back to rest of data
MZ_RT_filt_QC_avg <- cbind(MZ_RT_filt, avg_height_QC)

# check dimensions
dim(MZ_RT_filt_QC_avg)
```

    [1] 5406   74

Pull the name of your process blank, and make a new column that
indicates how many fold higher your peak height is in your average QC vs
your process blank.

``` r
# pull name of process blank 
grep("ProcessBlank", colnames(MZ_RT_filt), value = TRUE)
```

    [1] "ProcessBlank_pos_8"   "ProcessBlank2_pos_77"

``` r
# make a new column that has a value of how many fold higher peak height is in QCs as compared to PB
# then you can avg your PBs together and do the same thing
MZ_RT_filt_QC_avg$ProcessBlank_pos_8[is.na(MZ_RT_filt_QC_avg$ProcessBlank_pos_8)] <- 0
MZ_RT_filt_QC_avg$ProcessBlank2_pos_77[is.na(MZ_RT_filt_QC_avg$ProcessBlank2_pos_77)] <- 0

MZ_RT_filt_PB <- MZ_RT_filt_QC_avg %>% 
  mutate(avg_height_PB = (ProcessBlank2_pos_77 + ProcessBlank_pos_8)/2) %>%
  mutate(fold_higher_in_QC = avg_height_QC/avg_height_PB) %>%
  select(row_ID, mz_rt, mz, rt, avg_height_QC, avg_height_PB, fold_higher_in_QC)

head(MZ_RT_filt_PB)
```

      row_ID           mz_rt       mz        rt avg_height_QC avg_height_PB
    1    906 100.0394_0.4966 100.0394 0.4966192     12477.351      637.0004
    2    246 100.0748_0.3036 100.0748 0.3035967     13045.681      695.6903
    3   2700 100.0757_1.7369 100.0757 1.7369303     13219.717     5044.0909
    4   2770  100.112_1.8034 100.1120 1.8034459      7368.594        0.0000
    5   2472 100.1121_1.5856 100.1121 1.5856134      3181.897        0.0000
    6   2842 101.0072_1.8135 101.0072 1.8135151      6044.935     2296.0311
      fold_higher_in_QC
    1         19.587667
    2         18.752137
    3          2.620832
    4               Inf
    5               Inf
    6          2.632776

``` r
dim(MZ_RT_filt_PB)
```

    [1] 5406    7

We want to keep features that are at least 10x higher in QCs than
process blanks, and we also want to keep Infs, because an Inf indicates
that a feature absent in the process blanks (i.e., you get an Inf
because you’re trying to divide by zero).

``` r
# keep features that are present at least 10x higher in QCs vs PB
# or, keep NAs because those are absent in blank
PB_features_to_keep <- MZ_RT_filt_PB %>%
  filter(fold_higher_in_QC > 10 | is.infinite(fold_higher_in_QC)) 

dim(PB_features_to_keep)
```

    [1] 3289    7

How many features did we remove?

``` r
nrow(MZ_RT_filt_QC_avg) - nrow(PB_features_to_keep)
```

    [1] 2117

Wow removed a lot of garbage! This is great.

Bind back metdata.

``` r
MZ_RT_filt_PBremoved <- MZ_RT_filt_QC_avg %>%
  filter(mz_rt %in% PB_features_to_keep$mz_rt)
```

## Duplicates

Do we have any duplicate features?

``` r
get_dupes(MZ_RT_filt_PBremoved, mz_rt)
```

    No duplicate combinations found of: mz_rt

     [1] mz_rt                dupe_count           row_ID              
     [4] mz                   rt                   40TT2_pos_22        
     [7] FreshJuice5_pos_73   40TT4_pos_49         FreshJuice3_pos_64  
    [10] 70TT3_pos_14         FreshJuice0_pos_17   40TT6_pos_69        
    [13] 40TT_10M4_pos_26     FreshJuice4_pos_38   40TT1_pos_53        
    [16] 40TT5_pos_61         40TT3_pos_39         70TT2_pos_25        
    [19] 40TT_10M1_pos_18     FreshJuice6_pos_45   FreshJuice1_pos_24  
    [22] 40TT_10M3_pos_46     70TT7_pos_76         FreshJuice2_pos_54  
    [25] 70TT4_pos_55         40TT_10M2_pos_16     70TT5_pos_30        
    [28] 70TT6_pos_32         40TT_10M5_pos_40     40TT_10M6_pos_72    
    [31] 70TT1_pos_70         ProcessBlank_pos_8   ProcessBlank2_pos_77
    [34] QC6_pos_43           MEF1_pos_41          QC9_pos_67          
    [37] MEF_SS6_pos_75       MEF_SS5_pos_36       MEF4_pos_29         
    [40] MEF_SS4_pos_68       MEF_SS3_pos_52       MEF_SS1_pos_56      
    [43] MEF_SS2_pos_34       HPP1_pos_33          HPP4_pos_21         
    [46] HPP6_pos_31          MEF3_pos_44          MEF5_pos_50         
    [49] QC8_pos_59           MEF2_pos_60          HPP3_pos_27         
    [52] HPP2_pos_15          QC10_pos_74          QC7_pos_51          
    [55] HPP5_pos_71          QC2_pos_13           QC3_pos_20          
    [58] QC11_pos_78          QC5_pos_35           US1_pos_23          
    [61] UST2_pos_19          QC1_pos_12           UST1_pos_37         
    [64] US3_pos_47           QC4_pos_28           UST5_pos_63         
    [67] US5_pos_42           QC12_pos_79          US2_pos_48          
    [70] UST6_pos_62          US4_pos_65           UST3_pos_57         
    [73] US6_pos_66           UST4_pos_58          avg_height_QC       
    <0 rows> (or 0-length row.names)

In case there is duplicates, let’s remove them. The code below is not
currently running.

``` r
MZ_RT_filt_PBremoved_nodupes <- MZ_RT_filt_PBremoved %>%
  distinct(mz_rt, .keep_all = TRUE)

# how many duplicates did we remove?
nrow(MZ_RT_filt_PBremoved) - nrow(MZ_RT_filt_PBremoved_nodupes)

# check for dupes again, should be none
get_dupes(MZ_RT_filt_PBremoved_nodupes, mz_rt)

dim(MZ_RT_filt_PBremoved_nodupes)
```

Remove samples that we don’t need anymore.

``` r
MZ_RT_filt_PBremoved_extraremoved <- MZ_RT_filt_PBremoved %>%
  select(-ProcessBlank2_pos_77, -ProcessBlank_pos_8, -avg_height_QC, -FreshJuice0_pos_17)

dim(MZ_RT_filt_PBremoved_extraremoved)
```

    [1] 3289   70

Separate samples MEF x Pressure.

``` r
MZ_RT_filt_PBremoved_extraremoved_MEF <- MZ_RT_filt_PBremoved_extraremoved %>%
  select(-contains("HPP"), -contains("US"), -contains("UST"), -contains("70TT")) %>%
  select(-MEF_SS5_pos_36, -`40TT_10M1_pos_18`, -`40TT1_pos_53`, -`40TT2_pos_22` ,-`40TT3_pos_39`,-'40TT4_pos_49', -'40TT5_pos_61', -`40TT6_pos_69`)

dim(MZ_RT_filt_PBremoved_extraremoved_MEF)
```

    [1] 3289   37

``` r
glimpse(MZ_RT_filt_PBremoved_extraremoved_MEF)
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

``` r
MZ_RT_filt_PBremoved_extraremoved_pressure <- MZ_RT_filt_PBremoved_extraremoved %>%
  select(-contains("MEF"), -contains("MEF_SS"), -contains("40TT_10M")) %>% 
  select(-'70TT2_pos_25') 

dim(MZ_RT_filt_PBremoved_extraremoved_pressure)
```

    [1] 3289   52

``` r
glimpse(MZ_RT_filt_PBremoved_extraremoved_pressure)
```

    Rows: 3,289
    Columns: 52
    $ row_ID             <dbl> 906, 246, 2770, 2472, 2629, 750, 2370, 1895, 1299, …
    $ mz_rt              <chr> "100.0394_0.4966", "100.0748_0.3036", "100.112_1.80…
    $ mz                 <dbl> 100.0394, 100.0748, 100.1120, 100.1121, 101.0232, 1…
    $ rt                 <dbl> 0.4966192, 0.3035967, 1.8034459, 1.5856134, 1.67791…
    $ `40TT2_pos_22`     <dbl> NA, 9599.782, 9549.748, NA, 13851.662, 60896.620, 1…
    $ FreshJuice5_pos_73 <dbl> NA, 11359.062, 9083.935, NA, 5578.292, 53138.734, 5…
    $ `40TT4_pos_49`     <dbl> NA, 8453.625, 11118.646, NA, 15045.602, 70419.100, …
    $ FreshJuice3_pos_64 <dbl> NA, 6991.805, 10544.096, NA, 7092.719, 61444.938, 5…
    $ `70TT3_pos_14`     <dbl> NA, 10100.234, 8601.310, NA, 17391.098, 53557.120, …
    $ `40TT6_pos_69`     <dbl> NA, 8483.976, 7036.700, NA, 19088.390, 69080.830, 1…
    $ FreshJuice4_pos_38 <dbl> NA, 7958.649, 9853.188, 1311.417, 5236.963, 49287.3…
    $ `40TT1_pos_53`     <dbl> NA, 8169.814, 9138.023, NA, 14200.585, 67070.750, 7…
    $ `40TT5_pos_61`     <dbl> NA, 10176.031, 9042.211, NA, 20763.180, 65482.230, …
    $ `40TT3_pos_39`     <dbl> NA, 9242.485, 9607.523, NA, 17168.818, 64307.060, 9…
    $ FreshJuice6_pos_45 <dbl> NA, 9622.346, 9298.826, NA, 4703.250, 47034.375, 47…
    $ FreshJuice1_pos_24 <dbl> NA, 7768.952, 8173.136, NA, 3129.431, 51654.850, 47…
    $ `70TT7_pos_76`     <dbl> NA, 8460.284, 7806.664, NA, 19860.758, 55749.070, 1…
    $ FreshJuice2_pos_54 <dbl> NA, 8141.861, 6760.703, NA, 5957.557, 52588.240, 48…
    $ `70TT4_pos_55`     <dbl> NA, 10845.495, 7975.677, NA, 21994.710, 57688.992, …
    $ `70TT5_pos_30`     <dbl> NA, 7335.883, 5757.731, NA, 15414.226, 52490.035, 1…
    $ `70TT6_pos_32`     <dbl> NA, 7822.205, 7437.898, NA, 14236.880, 52129.023, 9…
    $ `70TT1_pos_70`     <dbl> NA, 8357.877, 7698.003, NA, 20410.920, 70099.880, 1…
    $ QC6_pos_43         <dbl> 12426.765, 16833.031, 7340.395, 3263.205, 18441.256…
    $ QC9_pos_67         <dbl> 13881.810, 14015.707, 9651.587, 4294.461, 20762.291…
    $ HPP1_pos_33        <dbl> NA, 8613.604, 8808.498, NA, 18136.605, 59198.688, 1…
    $ HPP4_pos_21        <dbl> NA, 10303.862, 12087.791, NA, 17066.062, 54223.258,…
    $ HPP6_pos_31        <dbl> NA, 7835.574, 8343.624, NA, 17386.156, 60233.426, 1…
    $ QC8_pos_59         <dbl> 10166.569, 11754.033, 7449.386, 3190.422, 16414.750…
    $ HPP3_pos_27        <dbl> NA, 8468.186, 10218.421, NA, 17978.572, 52360.305, …
    $ HPP2_pos_15        <dbl> NA, 9559.214, 8141.650, NA, 8241.522, 58958.440, 75…
    $ QC10_pos_74        <dbl> 12866.921, 14368.647, 8746.157, 3739.938, 21310.840…
    $ QC7_pos_51         <dbl> 14680.529, 16170.804, 11082.852, 1881.733, 15499.25…
    $ HPP5_pos_71        <dbl> NA, 8485.167, 7934.900, NA, 23159.232, 70392.320, 1…
    $ QC2_pos_13         <dbl> 12250.746, 14694.995, 5962.542, 3262.584, 15663.376…
    $ QC3_pos_20         <dbl> 13049.109, 17006.514, 6547.161, 3115.725, 17671.970…
    $ QC11_pos_78        <dbl> 14373.564, 14448.783, 9465.725, 4620.273, 19552.660…
    $ QC5_pos_35         <dbl> 16023.473, 12852.819, 6677.819, 3708.668, 16252.766…
    $ US1_pos_23         <dbl> NA, 8736.399, 5604.640, NA, 9204.489, 50559.918, 41…
    $ UST2_pos_19        <dbl> NA, 10806.581, 6584.571, NA, 14338.552, 53321.890, …
    $ QC1_pos_12         <dbl> 14941.828, 11572.707, 7141.102, 4035.967, 15206.507…
    $ UST1_pos_37        <dbl> NA, 10197.424, 5105.068, NA, 17689.795, 57432.855, …
    $ US3_pos_47         <dbl> NA, 10682.088, 7169.848, NA, 19364.693, 62903.734, …
    $ QC4_pos_28         <dbl> 12876.858, 13565.415, 8116.300, 3055.272, 15434.733…
    $ UST5_pos_63        <dbl> NA, 10110.027, 9069.041, NA, 25619.969, 63089.020, …
    $ US5_pos_42         <dbl> NA, 12408.161, 8666.140, NA, 18656.594, 51658.470, …
    $ QC12_pos_79        <dbl> 14667.275, 12310.267, 7610.508, 3196.210, 23105.748…
    $ US2_pos_48         <dbl> NA, 8647.708, 6084.618, NA, 17901.588, 70058.580, 9…
    $ UST6_pos_62        <dbl> NA, 8635.671, 10505.228, NA, 20837.086, 55077.312, …
    $ US4_pos_65         <dbl> NA, 8799.619, 9080.317, NA, 15373.005, 63397.580, 8…
    $ UST3_pos_57        <dbl> NA, 12310.286, 8217.092, NA, 19802.994, 58208.370, …
    $ US6_pos_66         <dbl> NA, 10080.784, 9462.883, 1174.036, 15593.913, 60810…
    $ UST4_pos_58        <dbl> NA, 8373.917, 7592.529, NA, 20930.805, 58288.344, 1…

# Save your file

Now you have a list of features present in your samples after filtering
for CV in QCs, and removing all the extraneous columns we added to help
us do this, along with removing any process blanks.

``` r
write_csv(MZ_RT_filt_PBremoved_extraremoved,
          "Post_filtering_pos_3289.csv")

write_csv(MZ_RT_filt_PBremoved_extraremoved_MEF,
          "Post_filtering_pos_3289_MEF.csv")

write_csv(MZ_RT_filt_PBremoved_extraremoved_pressure,
          "Post_filtering_pos_3289_pressure.csv")
```
