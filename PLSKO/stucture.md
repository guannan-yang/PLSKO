# Knockoff Steps and Functions Pipeline

## Genernal Input
- **Input: X**  
  Dimensions: `n x p`
  
- **Output: y**  
  Dimensions: `n x 1`
  
- **False Discovery Rate (FDR) Level**  
  Threshold: `q ∈ (0, 1]`


## Knockoff Framework Main Steps and corresponding functions included in this package

|                                            |                                  |                                              |                                                                                                    | **Pipeline Funs Provided^**   |                                    |
|--------------------------------------------|----------------------------------|----------------------------------------------|----------------------------------------------------------------------------------------------------|------------------------------|------------------------------------|
| **Mais Steps of Knockoff framework**       | **Function**                     | **Auxiliary Function**                       | **Output (Class)**                                                                                 | `plsko_filter()`, `plsAKO()` | `kofilter()`, `AKO_with_KO()`      |
| Step 1: Knockoff Variable Generation       | `plsko()`                        | `r_criterion()` <br>(for `ncomp` estimation) | A `n x p` matrix of knockoff variables                                                             | `Yes`                        | (Bring your own <br>knockoff vars) |
| Step 2: Importance Score Calculation (`W`) | (import from pkg `Knockoff`)     | -                                            | A vector of `p` or a `n_iter x p` matrix                                                           | `Yes`                        | `Yes`                              |
| Step 3: Knockoff Filtering                 | `KO_with_W()` <br>`AKO_with_W()` |                                              | A list (class "knockoff.result") with components: <br>`X`, `XKO`, `statistic`, `selected`, `index` | `Yes`                        | `Yes`                              |

### Additional Notes:
- ^The workflow includes generating knockoff variables, selecting models, and filtering.
- Graph results are S3 objects consistent with the knockoff method.
