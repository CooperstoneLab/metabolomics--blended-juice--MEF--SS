# Filtering - Untargeted Metabolomics Data - Negative Mode
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
metadata_juice <- read_csv(file = "Feature_list_all_data_neg.csv",
                      col_names = TRUE, # has headers
                      na = "0")
```

    Rows: 5119 Columns: 80
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    dbl (80): row ID, row m/z, row retention time, ConditioningQC_neg_5.mzML Pea...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# How many features do we have?
dim(metadata_juice)
```

    [1] 5119   80

``` r
# Look at beginning of the dataset

metadata_juice[1:8, 1:8]
```

    # A tibble: 8 × 8
      `row ID` `row m/z` `row retention time` ConditioningQC_neg_5.mzML Peak heigh…¹
         <dbl>     <dbl>                <dbl>                                  <dbl>
    1     1528      100.                0.898                                 12960.
    2      625      101.                0.338                                  4983.
    3      750      101.                0.337                                    NA 
    4      162      101.                0.291                                  1150.
    5      223      101.                0.294                                    NA 
    6      498      101.                0.321                                123496.
    7     1055      101.                0.478                                117850.
    8     1302      101.                0.740                                 60632.
    # ℹ abbreviated name: ¹​`ConditioningQC_neg_5.mzML Peak height`
    # ℹ 4 more variables: `40TT5_neg_61.mzML Peak height` <dbl>,
    #   `70TT1_neg_70.mzML Peak height` <dbl>,
    #   `40TT_10M1_neq_18.mzML Peak height` <dbl>,
    #   `ConditioningQC_neg_10.mzML Peak height` <dbl>

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

    [1] 0.300281 7.997691

``` r
# How many features do we have now?
dim(metadata_juice_RTfilt)
```

    [1] 4740   80

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
      row_ID mz_rt            mz    rt ConditioningQC_neg_5…¹ 40TT5_neg_61.mzML Pe…²
       <dbl> <chr>         <dbl> <dbl>                  <dbl>                  <dbl>
    1   1528 100.0402_0.8…  100. 0.898                 12960.                 12220.
    2    625 100.9439_0.3…  101. 0.338                  4983.                    NA 
    3    750 100.9512_0.3…  101. 0.337                    NA                     NA 
    4    498 101.0222_0.3…  101. 0.321                123496.                115426.
    5   1055 101.0238_0.4…  101. 0.478                117850.                115329.
    6   1302 101.0242_0.7…  101. 0.740                 60632.                 48460.
    7   2021 101.0243_1.7…  101. 1.73                  31236.                 27037.
    8   5612 101.0243_6.4…  101. 6.48                   6390.                  8019.
    # ℹ abbreviated names: ¹​`ConditioningQC_neg_5.mzML Peak height`,
    #   ²​`40TT5_neg_61.mzML Peak height`
    # ℹ 2 more variables: `70TT1_neg_70.mzML Peak height` <dbl>,
    #   `40TT_10M1_neq_18.mzML Peak height` <dbl>

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
     [4] "rt"                    "ConditioningQC_neg_5"  "40TT5_neg_61"         
     [7] "70TT1_neg_70"          "40TT_10M1_neq_18"      "ConditioningQC_neg_10"
    [10] "40TT_10M2_neq_16"      "ConditioningQC_neg_4"  "40TT3_neg_39"         
    [13] "40TT1_neg_53"          "70TT6_neg_32"          "FreshJuice1_neq_24"   
    [16] "70TT3_neq_14"          "ConditioningQC_neg_9"  "ConditioningQC_neg_3" 
    [19] "40TT2_neq_22"          "HPP1_neg_33"           "FreshJuice3_neg_64"   
    [22] "70TT4_neg_55"          "40TT_10M5_neg_40"      "40TT_10M4_neg_26"     
    [25] "FreshJuice4_neg_38"    "FreshJuice6_neg_45"    "FreshJuice5_neg_73"   
    [28] "40TT_10M6_neg_72"      "40TT4_neg_49"          "70TT7_neg_76"         
    [31] "70TT5_neg_30"          "70TT2_neq_25"          "40TT6_neg_69"         
    [34] "40TT_10M3_neg_46"      "FreshJuice_neq_17"     "SolventBlank4_neg_83" 
    [37] "SolventBlank2_neg_2"   "FreshJuice2_neg_54"    "ProcessBlank2_neg_77" 
    [40] "ProcessBlank1_neg_8"   "MEF_SS6_neg_75"        "SolventBlank1_neg_1"  
    [43] "MEF_SS4_neg_68"        "MEF_SS5_neg_36"        "MEF_SS1_neg_56"       
    [46] "MEF_SS2_neg_34"        "MEF_SS3_neg_52"        "MEF1_neg_41"          
    [49] "MEF2_neg_60"           "QC11_neg_78"           "QC4_neg_28"           
    [52] "QC10_neg_74"           "QC2_neg_13"            "MEF4_neg_29"          
    [55] "MEF5_neg_50"           "MEF3_neg_44"           "QC9_neg_67"           
    [58] "QC5_neg_35"            "QC6_neg_43"            "QC3_neq_20"           
    [61] "QC7_neg_51"            "QC8_neg_59"            "HPP4_neq_21"          
    [64] "US1_neq_23"            "QC12_neg_79"           "HPP6_neg_31"          
    [67] "HPP5_neg_71"           "HPP2_neq_15"           "HPP3_neg_27"          
    [70] "US2_neg_48"            "QC1_neg_12"            "US4_neg_65"           
    [73] "US3_neg_47"            "UST6_neg_62"           "UST1_neg_37"          
    [76] "UST2_neq_19"           "UST3_neg_57"           "US6_neg_66"           
    [79] "US5_neg_42"            "UST5_neg_63"           "UST4_neg_58"          

## Remove unwanted samples

Removing conditioning QCs and solvent blanks.

``` r
MZ_RT <- MZ_RT %>%
  select(-contains("Conditioning"), -contains("Solvent"))

colnames(MZ_RT)
```

     [1] "row_ID"               "mz_rt"                "mz"                  
     [4] "rt"                   "40TT5_neg_61"         "70TT1_neg_70"        
     [7] "40TT_10M1_neq_18"     "40TT_10M2_neq_16"     "40TT3_neg_39"        
    [10] "40TT1_neg_53"         "70TT6_neg_32"         "FreshJuice1_neq_24"  
    [13] "70TT3_neq_14"         "40TT2_neq_22"         "HPP1_neg_33"         
    [16] "FreshJuice3_neg_64"   "70TT4_neg_55"         "40TT_10M5_neg_40"    
    [19] "40TT_10M4_neg_26"     "FreshJuice4_neg_38"   "FreshJuice6_neg_45"  
    [22] "FreshJuice5_neg_73"   "40TT_10M6_neg_72"     "40TT4_neg_49"        
    [25] "70TT7_neg_76"         "70TT5_neg_30"         "70TT2_neq_25"        
    [28] "40TT6_neg_69"         "40TT_10M3_neg_46"     "FreshJuice_neq_17"   
    [31] "FreshJuice2_neg_54"   "ProcessBlank2_neg_77" "ProcessBlank1_neg_8" 
    [34] "MEF_SS6_neg_75"       "MEF_SS4_neg_68"       "MEF_SS5_neg_36"      
    [37] "MEF_SS1_neg_56"       "MEF_SS2_neg_34"       "MEF_SS3_neg_52"      
    [40] "MEF1_neg_41"          "MEF2_neg_60"          "QC11_neg_78"         
    [43] "QC4_neg_28"           "QC10_neg_74"          "QC2_neg_13"          
    [46] "MEF4_neg_29"          "MEF5_neg_50"          "MEF3_neg_44"         
    [49] "QC9_neg_67"           "QC5_neg_35"           "QC6_neg_43"          
    [52] "QC3_neq_20"           "QC7_neg_51"           "QC8_neg_59"          
    [55] "HPP4_neq_21"          "US1_neq_23"           "QC12_neg_79"         
    [58] "HPP6_neg_31"          "HPP5_neg_71"          "HPP2_neq_15"         
    [61] "HPP3_neg_27"          "US2_neg_48"           "QC1_neg_12"          
    [64] "US4_neg_65"           "US3_neg_47"           "UST6_neg_62"         
    [67] "UST1_neg_37"          "UST2_neq_19"          "UST3_neg_57"         
    [70] "US6_neg_66"           "US5_neg_42"           "UST5_neg_63"         
    [73] "UST4_neg_58"         

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

    [1] 4740   73

``` r
MZ_RT_QCs <- MZ_RT %>%
  select(mz_rt, contains("QC")) %>% # select QCs
  filter(rowSums(is.na(.)) <= 1) # remove rows that have 1 or more NAs
```

``` r
# check dimensions of QCs filtered df
dim(MZ_RT_QCs)
```

    [1] 4654   13

``` r
# how many features got removed with this filtering?
nrow(MZ_RT) - nrow(MZ_RT_QCs)
```

    [1] 86

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

    [1] 201

## Merge back the rest of the data

MZ_RT_QCs_CVfilt only contains the QCs, We want to keep only the rows
that are present in this df, and then merge back all of the other
samples present in MZ_RT. We will do this by creating a vector that has
the mz_rt features we want to keep, and then using `filter()` and `%in%`
to keep only features that are a part of this list.

``` r
dim(MZ_RT_QCs_CVfilt)
```

    [1] 4453   14

``` r
dim(MZ_RT)
```

    [1] 4740   73

``` r
# make a character vector of the mz_rt features we want to keep
# i.e., the ones that passed our previous filtering steps

features_to_keep <- as.character(MZ_RT_QCs_CVfilt$mz_rt)

MZ_RT_filt <- MZ_RT %>%
  filter(mz_rt %in% features_to_keep)

dim(MZ_RT_filt)
```

    [1] 4453   73

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

    [1] "ProcessBlank2_neg_77" "ProcessBlank1_neg_8" 

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

    [1] 4453   74

Pull the name of your process blank, and make a new column that
indicates how many fold higher your peak height is in your average QC vs
your process blank.

``` r
# pull name of process blank 
grep("ProcessBlank", colnames(MZ_RT_filt), value = TRUE)
```

    [1] "ProcessBlank2_neg_77" "ProcessBlank1_neg_8" 

``` r
# make a new column that has a value of how many fold higher peak height is in QCs as compared to PB
# then you can avg your PBs together and do the same thing
MZ_RT_filt_QC_avg$ProcessBlank1_neg_8[is.na(MZ_RT_filt_QC_avg$ProcessBlank1_neg_8)] <- 0
MZ_RT_filt_QC_avg$ProcessBlank2_neg_77[is.na(MZ_RT_filt_QC_avg$ProcessBlank2_neg_77)] <- 0

MZ_RT_filt_PB <- MZ_RT_filt_QC_avg %>% 
  mutate(avg_height_PB = (ProcessBlank2_neg_77 + ProcessBlank1_neg_8)/2) %>%
  mutate(fold_higher_in_QC = avg_height_QC/avg_height_PB) %>%
  select(row_ID, mz_rt, mz, rt, avg_height_QC, avg_height_PB, fold_higher_in_QC)

head(MZ_RT_filt_PB)
```

      row_ID           mz_rt       mz        rt avg_height_QC avg_height_PB
    1   1528 100.0402_0.8975 100.0402 0.8975327     11200.859         0.000
    2    625 100.9439_0.3382 100.9439 0.3381544      5282.121         0.000
    3    498 101.0222_0.3207 101.0222 0.3207379    107932.718         0.000
    4   1055 101.0238_0.4781 101.0238 0.4781379    109764.404      6447.208
    5   1302 101.0242_0.7398 101.0242 0.7397535     47894.475      6922.726
    6   2021 101.0243_1.7271 101.0243 1.7270532     30130.920      7121.128
      fold_higher_in_QC
    1               Inf
    2               Inf
    3               Inf
    4         17.025107
    5          6.918441
    6          4.231200

``` r
dim(MZ_RT_filt_PB)
```

    [1] 4453    7

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

    [1] 2997    7

How many features did we remove?

``` r
nrow(MZ_RT_filt_QC_avg) - nrow(PB_features_to_keep)
```

    [1] 1456

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
     [4] mz                   rt                   40TT5_neg_61        
     [7] 70TT1_neg_70         40TT_10M1_neq_18     40TT_10M2_neq_16    
    [10] 40TT3_neg_39         40TT1_neg_53         70TT6_neg_32        
    [13] FreshJuice1_neq_24   70TT3_neq_14         40TT2_neq_22        
    [16] HPP1_neg_33          FreshJuice3_neg_64   70TT4_neg_55        
    [19] 40TT_10M5_neg_40     40TT_10M4_neg_26     FreshJuice4_neg_38  
    [22] FreshJuice6_neg_45   FreshJuice5_neg_73   40TT_10M6_neg_72    
    [25] 40TT4_neg_49         70TT7_neg_76         70TT5_neg_30        
    [28] 70TT2_neq_25         40TT6_neg_69         40TT_10M3_neg_46    
    [31] FreshJuice_neq_17    FreshJuice2_neg_54   ProcessBlank2_neg_77
    [34] ProcessBlank1_neg_8  MEF_SS6_neg_75       MEF_SS4_neg_68      
    [37] MEF_SS5_neg_36       MEF_SS1_neg_56       MEF_SS2_neg_34      
    [40] MEF_SS3_neg_52       MEF1_neg_41          MEF2_neg_60         
    [43] QC11_neg_78          QC4_neg_28           QC10_neg_74         
    [46] QC2_neg_13           MEF4_neg_29          MEF5_neg_50         
    [49] MEF3_neg_44          QC9_neg_67           QC5_neg_35          
    [52] QC6_neg_43           QC3_neq_20           QC7_neg_51          
    [55] QC8_neg_59           HPP4_neq_21          US1_neq_23          
    [58] QC12_neg_79          HPP6_neg_31          HPP5_neg_71         
    [61] HPP2_neq_15          HPP3_neg_27          US2_neg_48          
    [64] QC1_neg_12           US4_neg_65           US3_neg_47          
    [67] UST6_neg_62          UST1_neg_37          UST2_neq_19         
    [70] UST3_neg_57          US6_neg_66           US5_neg_42          
    [73] UST5_neg_63          UST4_neg_58          avg_height_QC       
    <0 rows> (or 0-length row.names)

No duplicates found. In case there is duplicates, let’s remove them. The
code below is not currently running.

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
  select(-ProcessBlank2_neg_77, -ProcessBlank1_neg_8, -avg_height_QC, -FreshJuice_neq_17)

dim(MZ_RT_filt_PBremoved_extraremoved)
```

    [1] 2997   70

Separate samples MEF x Pressure.

``` r
MZ_RT_filt_PBremoved_extraremoved_MEF <- MZ_RT_filt_PBremoved_extraremoved %>%
  select(-contains("HPP"), -contains("US"), -contains("UST"), -contains("70TT")) %>%
  select(-MEF_SS5_neg_36, -`40TT_10M1_neq_18`, -`40TT1_neg_53`, -`40TT2_neq_22` ,-`40TT3_neg_39`,-'40TT4_neg_49', -'40TT5_neg_61', -`40TT6_neg_69`)

dim(MZ_RT_filt_PBremoved_extraremoved_MEF)
```

    [1] 2997   37

``` r
glimpse(MZ_RT_filt_PBremoved_extraremoved_MEF)
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

``` r
MZ_RT_filt_PBremoved_extraremoved_pressure <- MZ_RT_filt_PBremoved_extraremoved %>%
  select(-contains("MEF"), -contains("MEF_SS"), -contains("40TT_10M")) %>% 
  select(-'70TT2_neq_25') 

dim(MZ_RT_filt_PBremoved_extraremoved_pressure)
```

    [1] 2997   52

``` r
glimpse(MZ_RT_filt_PBremoved_extraremoved_pressure)
```

    Rows: 2,997
    Columns: 52
    $ row_ID             <dbl> 1528, 625, 498, 1055, 345, 1117, 1002, 630, 854, 12…
    $ mz_rt              <chr> "100.0402_0.8975", "100.9439_0.3382", "101.0222_0.3…
    $ mz                 <dbl> 100.0402, 100.9439, 101.0222, 101.0238, 101.0294, 1…
    $ rt                 <dbl> 0.8975327, 0.3381544, 0.3207379, 0.4781379, 0.32149…
    $ `40TT5_neg_61`     <dbl> 12219.674, NA, 115425.840, 115329.305, 7823.426, 84…
    $ `70TT1_neg_70`     <dbl> 14698.048, NA, 91239.960, 119028.520, 7469.719, 824…
    $ `40TT3_neg_39`     <dbl> 12310.439, NA, 87698.370, 122651.580, 6439.374, 983…
    $ `40TT1_neg_53`     <dbl> 10453.251, NA, 117572.780, 113535.305, 9150.244, 93…
    $ `70TT6_neg_32`     <dbl> 10093.323, NA, 110411.610, 112835.530, 9580.913, 93…
    $ FreshJuice1_neq_24 <dbl> 12023.159, NA, 88496.350, 112586.830, 6805.068, 854…
    $ `70TT3_neq_14`     <dbl> 11960.279, NA, 114490.620, 128249.670, 8802.522, 10…
    $ `40TT2_neq_22`     <dbl> 13026.044, NA, 99735.336, 112527.290, 7357.439, 771…
    $ HPP1_neg_33        <dbl> 16526.880, NA, 86278.945, 120722.195, 5406.450, 938…
    $ FreshJuice3_neg_64 <dbl> 13652.454, NA, 96339.530, 112846.810, 7138.783, 865…
    $ `70TT4_neg_55`     <dbl> 13600.727, NA, 91928.414, 116925.120, 6745.480, 836…
    $ FreshJuice4_neg_38 <dbl> 13704.267, NA, NA, 115320.030, NA, 8136.055, 6654.9…
    $ FreshJuice6_neg_45 <dbl> 11589.513, NA, 98065.090, 111706.766, 7787.599, 822…
    $ FreshJuice5_neg_73 <dbl> 13904.639, NA, 96784.640, 108390.480, 7504.919, 759…
    $ `40TT4_neg_49`     <dbl> 11945.086, NA, 100951.570, 114384.470, 8183.823, 84…
    $ `70TT7_neg_76`     <dbl> 12456.034, NA, 95837.830, 115806.480, 7758.618, 808…
    $ `70TT5_neg_30`     <dbl> 13512.075, NA, 88549.900, 124891.860, 6005.645, 939…
    $ `40TT6_neg_69`     <dbl> 11734.688, NA, 89414.914, 113368.590, 6610.147, 830…
    $ FreshJuice2_neg_54 <dbl> 13125.755, NA, 87831.695, 105103.160, 7000.599, 743…
    $ QC11_neg_78        <dbl> 13470.114, 5929.262, 120178.180, 125090.984, 8608.3…
    $ QC4_neg_28         <dbl> 13730.859, 5456.804, 108522.360, 118438.914, 7660.1…
    $ QC10_neg_74        <dbl> 12069.027, 4993.207, 109453.290, 118015.520, 7795.1…
    $ QC2_neg_13         <dbl> 13325.600, 6165.458, 119287.766, 125860.680, 9659.3…
    $ QC9_neg_67         <dbl> 11914.959, 6696.042, 92050.380, 113209.550, 6111.62…
    $ QC5_neg_35         <dbl> 11789.217, 7697.496, 101195.160, 112213.750, 6851.1…
    $ QC6_neg_43         <dbl> 11772.166, 4621.732, 176979.610, 116165.266, 11990.…
    $ QC3_neq_20         <dbl> 11308.223, 3726.545, 119111.340, 121476.016, 9253.6…
    $ QC7_neg_51         <dbl> 11150.521, 3754.003, 109563.060, 115072.195, 7161.4…
    $ QC8_neg_59         <dbl> 11563.439, 6900.274, 114135.940, 119587.280, 8842.6…
    $ HPP4_neq_21        <dbl> 13258.536, NA, 109516.766, 114320.500, 8561.948, 89…
    $ US1_neq_23         <dbl> 9220.965, NA, 120037.930, 96670.516, 8141.219, 7465…
    $ QC12_neg_79        <dbl> 12009.428, 5218.626, 112950.984, 112701.120, 7004.2…
    $ HPP6_neg_31        <dbl> 11201.091, 1528.570, 110677.625, 118004.700, 9369.7…
    $ HPP5_neg_71        <dbl> 13716.359, NA, 107396.890, 121256.984, NA, 7408.554…
    $ HPP2_neq_15        <dbl> 13398.194, NA, 114677.550, 121721.320, 8692.591, 96…
    $ HPP3_neg_27        <dbl> 16690.006, NA, 107935.720, 119360.470, 8404.860, 92…
    $ US2_neg_48         <dbl> 9830.868, NA, 100541.055, 109872.266, 7212.390, 802…
    $ QC1_neg_12         <dbl> 11507.536, 7507.896, 119697.086, 129105.930, 8513.6…
    $ US4_neg_65         <dbl> 11588.510, NA, 92703.680, 112420.260, 7097.869, 950…
    $ US3_neg_47         <dbl> 11859.217, NA, 88344.805, 115693.160, 6168.087, 940…
    $ UST6_neg_62        <dbl> 14063.892, NA, 89911.234, 113612.400, 5847.746, 838…
    $ UST1_neg_37        <dbl> 12078.936, NA, 98425.380, 108826.970, 8379.202, 776…
    $ UST2_neq_19        <dbl> 9850.374, NA, 106590.840, 107984.380, NA, 8615.376,…
    $ UST3_neg_57        <dbl> 10960.633, NA, 102329.336, 112129.560, 6566.236, 74…
    $ US6_neg_66         <dbl> 10767.037, NA, 96531.305, 114152.660, 7065.190, 740…
    $ US5_neg_42         <dbl> 14527.405, NA, 81132.690, 114641.630, 6365.109, 884…
    $ UST5_neg_63        <dbl> 13444.101, NA, 99419.030, 115687.780, 8132.580, 821…
    $ UST4_neg_58        <dbl> 12354.296, 1172.319, 82394.960, 114271.240, 5364.44…

# Save your file

Now you have a list of features present in your samples after filtering
for CV in QCs, and removing all the extraneous columns we added to help
us do this, along with removing any process blanks.

``` r
write_csv(MZ_RT_filt_PBremoved_extraremoved,
          "Post_filtering_neg_2997.csv")

write_csv(MZ_RT_filt_PBremoved_extraremoved_MEF,
          "Post_filtering_neg_2997_MEF.csv")

write_csv(MZ_RT_filt_PBremoved_extraremoved_pressure,
          "Post_filtering_neg_2997_pressure.csv")
```
