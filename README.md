# PLSKO
Package and source codes used in paper: `PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection'

# Knockoff Steps and Functions Pipeline

## Genernal Input
- **Input: predictor matrix X**  
  Dimensions: `n x p`
  
- **Input: response vector y**  
  Dimensions: `n x 1`
  
- **Input: Target False Discovery Rate (FDR) Level**  
  Threshold: `q ∈ (0, 1]`

- **Output: `knockoff.result` object with Selected variables**


## Knockoff framework main steps and corresponding functions included in this package

|                                                |                                  |                                              |                                                                                                     | **Pipeline Funs Provided^**  |                                    |
|------------------------------------------------|----------------------------------|----------------------------------------------|-----------------------------------------------------------------------------------------------------|------------------------------|------------------------------------|
| **Mais Steps of Knockoff framework**           | **Function**                     | **Auxiliary Function**                       | **Output (_Class_)**                                                                                | `plsko_filter()`, `plsAKO()` | `ko_filter()`, `AKO_with_KO()`     |
| **Step 1: Knockoff Variable Generation**       | `plsko()`                        | `r_criterion()` <br>(for `ncomp` estimation) | A `n x p` _matrix_ of knockoff variables                                                            | :heavy_check_mark:           | (Bring your own <br>knockoff vars) |
| **Step 2: Importance Score Calculation (`W`)** | (import from pkg `Knockoff`)     | -                                            | A _vector_ of `p` or a `n_ko x p` _matrix_                                                        | :heavy_check_mark:           | :heavy_check_mark:                 |
| **Step 3: Knockoff Filtering and Variable Selection**                 | `KO_with_W()` <br>`AKO_with_W()` | -                                             | A list (class _"knockoff.result"_ or _"AKO.result"_) with components: <br> `statistic`, `selected`, `ako.selected` | :heavy_check_mark:           | :heavy_check_mark:                 |

### Additional Notes:
- ^The workflow includes generating knockoff variables, importance score calculation, and filtering.

## Vignette 
[View the vignette](https://github.com/guannan-yang/PLSKO/blob/main/docs/PLSKO.html)
