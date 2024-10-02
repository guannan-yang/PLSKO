
# PLSKO
Package and source codes used in paper, 'PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection'(https://www.biorxiv.org/content/10.1101/2024.08.06.606935v1)

<img src="https://github.com/guannan-yang/PLSKO/blob/main/docs/plsko.png" width=25% height=25%>

# Installation
You can install the development version of the package from GitHub using the following code:
```{r installation}
# install.packages("devtools")
devtools::install_github("guannan-yang/PLSKO/PLSKO", quiet = TRUE)
```
If warnings of installing dependencies appears, please install the dependencies manually by running the following code:
```{r install dependencies}
# install.packages(c("knockoff","progress", "parallel", "doParallel", "foreach"))
#
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("mixOmics")
```

# Vignette 
[View the vignette](https://guannan-yang.github.io/PLSKO/PLSKO.html)

# Knockoff Steps and Functions Pipeline

## Basic Input
- **Input: predictor matrix X**  
  Dimensions: `n x p`
  
- **Input: response vector y**  
  Length of `n`
  
- **Input: Target False Discovery Rate (FDR) Level**  
  Threshold: `q ∈ (0, 1]`

- **Output: `knockoff.result` object with Selected variables**


## Main functions included in this package

|                                                |                                  |                                              |                                                                                                     | **Pipeline Functions Provided^**  |                                    |
|------------------------------------------------|----------------------------------|----------------------------------------------|-----------------------------------------------------------------------------------------------------|------------------------------|------------------------------------|
| **Main Steps of Knockoff framework**           | **Function**                     | **Auxiliary Function**                       | **Output (_Class_)**                                                                                | `plsko_filter()`, `plsAKO()` | `ko_filter()`, `AKO_with_KO()`     |
| **Step 1: Knockoff Variable Generation**       | `plsko()`                        | `r_criterion()` <br>(for `ncomp` estimation) | A `n x p` _matrix_ of knockoff variables                                                            | :heavy_check_mark:           | (Bring your own <br>knockoff vars) |
| **Step 2: Importance Score Calculation (`W`)** | (import from pkg `Knockoff`)     | -                                            | A _vector_ of `p` or a `n_ko x p` _matrix_                                                        | :heavy_check_mark:           | :heavy_check_mark:                 |
| **Step 3: Knockoff Filtering and Variable Selection**                 | `KO_with_W()` <br>`AKO_with_W()` | -                                             | A list (class _"knockoff.result"_ or _"AKO.result"_) with components: <br> `statistic`, `selected`, `ako.selected` | :heavy_check_mark:           | :heavy_check_mark:                 |

[View the Knockoff Framework Diagram](https://github.com/guannan-yang/PLSKO/blob/main/docs/diagram-ko.pdf)

### Additional Notes:
- ^The workflow includes generating knockoff variables, importance score calculation, and filtering.


