# PLSKO
Package (in preparation) and source codes used in paper: `PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection'

# Knockoff Steps and Functions Pipeline

## Genernal Input
- **Input: X**  
  Dimensions: `n x p`
  
- **Output: y**  
  Dimensions: `n x 1`
  
- **False Discovery Rate (FDR) Level**  
  Threshold: `q ∈ (0, 1]`


## Knockoff Framework Main Steps and corresponding functions included in this package

|                                                |                                  |                                              |                                                                                                   | **Pipeline Funs Provided^**  |                                    |
|------------------------------------------------|----------------------------------|----------------------------------------------|---------------------------------------------------------------------------------------------------|------------------------------|------------------------------------|
| **Mais Steps of Knockoff framework**           | **Function**                     | **Auxiliary Function**                       | **Output (Class)**                                                                                | `plsko_filter()`, `plsAKO()` | `ko_filter()`, `AKO_with_KO()`     |
| **Step 1: Knockoff Variable Generation**       | `plsko()`                        | `r_criterion()` <br>(for `ncomp` estimation) | A `n x p` matrix of knockoff variables                                                            | - [x]                        | (Bring your own <br>knockoff vars) |
| **Step 2: Importance Score Calculation (`W`)** | (import from pkg `Knockoff`)     | -                                            | A vector of `p` or a `n_iter x p` matrix                                                          | - [x]                        | - [x]                              |
| **Step 3: Knockoff Filtering**                 | `KO_with_W()` <br>`AKO_with_W()` |                                              | A list (class "knockoff.result") with components: <br>`X`, `Xk`, `statistic`, `selected`, `index` | - [x]                        | - [x]                              |

### Additional Notes:
- ^The workflow includes generating knockoff variables, selecting models, and filtering.
