# Probabilistic prediction with vine copulas

Code of the Medium post *"Your Model Says 0.73. But How Sure Is It?"*: a breast-cancer classifier and a
California house-price model whose outputs are calibrated probabilities and full predictive distributions,
both built on one fortvine vine copula.

| file | what it does |
|---|---|
| `vinepred.py` | the models: smooth KDE margins, `VineClassifier` (the class label as a discrete vine variable), `VineRegressor` (predictive density on a grid via `Vinecop.pdf_grid`); `python vinepred.py` runs a self-check |
| `breast_cancer.py` | 10-fold cross-validation on the Wisconsin diagnostic data against naive Bayes, logistic regression and gradient boosting; 14 metrics |
| `california.py` | 80/20 split on California housing against quantile boosting and a Gaussian linear model; NSE, KGE, CRPS, coverage and more |
| `figures.py` | the seven figures of the post, into `img/` |

```bash
pip install fortvine scikit-learn matplotlib scipy
python breast_cancer.py      # writes bc_results.json, bc_pred.npz
python california.py         # downloads the data once; writes ca_results.json, ca_pred.npz
python figures.py
```

The scripts use up to 12 threads (`num_threads=12`); lower it on a smaller machine.
