# ActivityProfileR
---

## Overview
---

ActivityProfileR is an implementation of the functional ANOVA test for activity data described in *Nonparametric comparisons of activity profiles from wearable device data* by Hsin-wen Chang. A presentation on the original paper can be found [here](https://www.ima.umn.edu/materials/2019-2020/DW9.16-17.19/28237/talk_Minneapolis.pdf). Once your data is in the correct format (described below), you can pass it into `EL_test()`, and it will conduct the statistical test for you. 

ActivityProfileR uses the tidyverse libraries for data cleaning and manipulation. There is also a paralellized version of the test in `EL_test_parallel()` available using multidplyr.

The statistical test can take a while to run, so parallelizaiton is recommended.

## Installation
---

```
# install.packages("devtools")
devtools::install_github("ActivityProfileR")
```

## Usage
---

To use ActivityProfileR, you must first get your activity data into the correct format. The `EL_test()` function expects the activity data for each subject to be nested, so that each subject will only have one row. 

You can use `prepare_sim_data()` function to create an example of the data `EL_test()` expects.

```
library(tidyverse)
library(rootSolve)
library(ActivityProfileR)

# Generate data for 3 groups with different characteristics
params_by_group = list(
  list(n_subj = 70, n_points = 1000, 
       scaling = 300, thetas = c(-0.8, 2.0, 2.31)),
  list(n_subj = 100, n_points = 1000, 
       scaling = 300, thetas = c(-0.402, 1.0, 1.65)),
  list(n_subj = 130, n_points = 1000, 
       scaling = 300, thetas = c(-0.201, 0.5, 1.18))
)

sim_data = prepare_sim_data(params_by_group)
```

```
head(sim_data)
## # A tibble: 6 x 3
##   subject_id group Xt           
##        <int> <dbl> <list>       
## 1          1     1 <dbl [1,000]>
## 2          2     1 <dbl [1,000]>
## 3          3     1 <dbl [1,000]>
## 4          4     1 <dbl [1,000]>
## 5          5     1 <dbl [1,000]>
## 6          6     1 <dbl [1,000]>
```

Now that the data is created, it can be passed into the `EL_test()` function to perform the test.

```
# Conduct the functional ANOVA test
el = EL_test(sim_data)
```

The test uses bootstrap and uses optimization, so it can take a while to run, about 4 minutes. We can also use the parallelized version of the test `EL_test_parallel()`. This version of the function expects a number of cores that it should use.

```
# Determine the number of cores to use in parallelized version of functions
n_cores = parallel::detectCores() - 2
par_el = EL_test_parallel(sim_data, n_cores = n_cores)
```

Parallelization can cut down calculation time by about half. See vignettes("demo") for more details.