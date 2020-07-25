# ActivityProfileR

## Overview

ActivityProfileR is an implementation of the functional ANOVA test for activity data described in *Nonparametric comparisons of activity profiles from wearable device data* by Hsin-wen Chang and Professor Ian McKeague. A presentation on the original paper can be found [here](https://www.ima.umn.edu/materials/2019-2020/DW9.16-17.19/28237/talk_Minneapolis.pdf).  

## Installation

```
# install.packages("remotes")
remotes::install_github("thecbp/ActivityProfileR")
```

## Usage

To use `ActivityProfileR`, the package expects activity data in a particular format. You can see an example of this data format using the `prepare_sim_data()` function. The `prepare_sim_data()` expects a list of lists, where each inner lists contains parameters for a group you want to appear in the data. Activity data is generated through an Ornstein-Uhlenbeck process (implemented in the `ESGtoolkit` package).

```
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

data = prepare_sim_data(params_by_group)
```

```
head(data)
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

Your data should have the following information:

- `subject_id` identifies a subject's ID within their group. For example, the number `1` will identify the first person in a group, so there should be multiple `1` in this column, and so on.
- `group` identifies the group ID that a subject is in
- `Xt` contains the activity data for each subject, nested for each 

After getting your data into this format, you will need to use the `add_helper_columns()` and `add_activity_profile()` functions to further process the data. `add_helper_columns()` helps add some useful columns for later analysis, and `add_activity_profile()` creates a column that converts the raw activity data into the activity profile that is the crux of the test.

```
processed_data = data %>% 
    add_helper_columns(group_col = group) %>% 
    add_activity_profile(group_col = group, 
                         activity_col = Xt)
```

You will need to provide the names of the columns that contain this corresponding information in your data (subject ID, group ID, activity data). Once your data is in this format, it can be passed into the `EL_test()` function to perform the test.

```
# Test uses bootstrapping, so set seed for consistency
set.seed(1)

# Conduct the functional ANOVA test
el = processed_data %>% 
  EL_test(id_col = subject_id,
          group_col = group,
          activity_col = Xt)
```

The test uses bootstrap and uses optimization, so it can take a while to run. See `vignettes("demo")` for more details.