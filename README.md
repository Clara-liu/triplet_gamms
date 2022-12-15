# GAMM analysis for finding onsets in minimal triplets
The script here contains the function used for the GAMM analysis in [Liu, Xu and Hsieh (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0095447021000917).
## Installation
1. Clone the current repo. Alternatively if you do not use git, download the `function.R` script.
2. Install the required R packages in your R terminal. The scripted is tested with version 2.4.1 for `itsadug` and 1.8-39 for `mgcv`.
```
install.packages(c('mgcv', 'itsadug'))
```
3. Run the content in `function.R`.
## Preparation
### Stimuli list
Create a stimuli list in R. An example is given here in `stimuli_example.R`. Note that `t1` represents triplet 1. The triplet key should correspond to what is used in your dataset.

The first two words should be the vowel minimal pair and the second the consonant minimal pair.
### Data preparation
The data need to be arranged in an R list with each triplet being a list element. e.g. 6 R dataframes as the 6 element in a large list. Each dataframe should consist of data for one triplet.

An example for one dataframe is given here in `data_example.csv`. The dataframes should contain the following columns:
- *Time*: Time stamp for each row.
- *Triplet*: Information indicating which triplet the current time series belong to. This need to contain the same labels as what is used as the keys in the stimuli list.
- *Rep*: Which repetition the current time series belong to. This is used as a random effect in the GAMMs.
- *Word*: Which word the current time series belong to. This need to contain the same word labels as the values in the stimuli list.
- Your measurement column such as `F2`, `ULx` etc.
- *Speaker*: Which speaker the current time series belong to.

Because GAMM requires for all the time series to contain the same number of samples, each time series should be time normalised or trimmed to have the same number of rows. The example data here contains 96 time points per word sequence.
## Usage
4 arguments are needed to run the function:
- *df_list*: The list of dataframes mentioned above.
- *k_*: The k parameter used in the GAMMs. You can determine this by manually modelling for one triplet and checking with the `gam.check()` function in the `mgcv` package.
- *min_sig_length*: How long the divergence need to be significant for the onset to be recorded. 40 ms was used in [Liu et al (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0095447021000917).
- *measurement*: Which measurement to use. e.g. "F2".
### Example
```R
measurement<- 'F2'
sig_interval<- 0.04
k<-15
source('stimuli_example.R')
source('function.R')
# overview of the data list
str(s1_dfs, max.level = 1)
```
```R
List of 6
 $ :'data.frame':	2790 obs. of  32 variables:
 $ :'data.frame':	2790 obs. of  32 variables:
 $ :'data.frame':	2790 obs. of  32 variables:
 $ :'data.frame':	2790 obs. of  32 variables:
 $ :'data.frame':	2790 obs. of  32 variables:
 $ :'data.frame':	2790 obs. of  32 variables:
```
```R
# overview of one dataframe in the data list
str(s1_dfs[[1]])
```
```R
'data.frame':	2790 obs. of  32 variables:
 $ File       : Factor w/ 30 levels "lailiwei01","lailiwei02",..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Interval   : Factor w/ 4 levels "lai","li","lu",..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Time       : num  0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 ...
 $ ULx        : num  12.1 12.1 12 12 12 ...
 $ ULy        : num  -0.644 -0.628 -0.612 -0.596 -0.564 ...
 $ ULz        : num  -0.1737 -0.1425 -0.1113 -0.0806 -0.0899 ...
 $ LLx        : num  7.67 7.6 7.53 7.47 7.4 ...
 $ LLy        : num  -2.54 -2.48 -2.41 -2.35 -2.3 ...
 $ LLz        : num  -27.9 -28 -28.1 -28.2 -28.2 ...
 $ TTx        : num  -12.6 -12.5 -12.5 -12.4 -12.3 ...
 $ TTy        : num  -2.47 -2.45 -2.43 -2.41 -2.33 ...
 $ TTz        : num  0.936 0.931 0.926 0.921 0.915 ...
 $ TBx        : num  -31.7 -31.7 -31.6 -31.6 -31.6 ...
 $ TBy        : num  -3.68 -3.68 -3.68 -3.68 -3.72 ...
 $ TBz        : num  -4.36 -4.63 -4.89 -5.16 -5.42 ...
 $ TDx        : num  -44.4 -44.4 -44.4 -44.3 -44.3 ...
 $ TDy        : num  -0.899 -0.905 -0.911 -0.917 -0.93 ...
 $ TDz        : num  -5.65 -5.97 -6.3 -6.62 -6.92 ...
 $ JAWx       : num  -1.79 -1.78 -1.78 -1.77 -1.77 ...
 $ JAWy       : num  -1.72 -1.71 -1.7 -1.69 -1.68 ...
 $ JAWz       : num  -25.9 -26 -26 -26 -26.1 ...
 $ F1         : num  350 351 353 354 356 ...
 $ F2         : num  0.318 0.289 0.261 0.233 0.238 ...
 $ F3         : num  3044 3093 3143 3193 3212 ...
 $ F2_F3      : num  2375 2397 2418 2440 2468 ...
 $ Speaker    : Factor w/ 1 level "s1": 1 1 1 1 1 1 1 1 1 1 ...
 $ LP         : num  -0.516 -0.546 -0.575 -0.604 -0.635 ...
 $ Rep        : Factor w/ 10 levels "01","02","03",..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Word       : Factor w/ 3 levels "lailiwei","lailuwei",..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Int        : Factor w/ 2 levels "a","b": 1 1 1 1 1 1 1 1 1 1 ...
 $ start.event: logi  TRUE FALSE FALSE FALSE FALSE FALSE ...
 $ Triplet    : chr  "t1" "t1" "t1" "t1" ...
```
```R
# run the analysis
s1_f2_diverge<- getDivergeData(s1_dfs, k, stimuli_list, measurement, sig_interval)
# overview of the result
str(s1_f2_diverge)
```
```R
'data.frame':	12 obs. of  5 variables:
 $ Type      : chr  "v" "c" "v" "c" ...
 $ Diverge   : num  0.232 0.218 0.204 0.269 0.269 ...
 $ Triplet   : chr  "1" "1" "2" "2" ...
 $ Syl_struct: chr  "t" "t" "t" "t" ...
 $ Speaker   : chr  "s1" "s1" "s1" "s1" ...
```
