# Data dictionary — Breast Cancer Survival

Source: the Kaggle public dataset *"Breast Cancer Survival"* (originally an
SPSS `.sav` file). 1,207 breast-cancer patients, 9 columns. The original
`.sav` is kept here for provenance; `breast_cancer_survival.csv` is a clean
CSV export used as the primary runtime input.

| Variable  | Type        | Description                                              | Missing |
|-----------|-------------|----------------------------------------------------------|---------|
| `id`      | int         | Patient identifier (dropped before modelling)            | 0       |
| `age`     | numeric     | Age at diagnosis (years)                                 | 0       |
| `pathsize`| numeric     | Pathologic tumor size (cm)                               | 86      |
| `lnpos`   | numeric     | Number of positive axillary lymph nodes                  | 0       |
| `histgrad`| ordinal     | Histologic grade: 1 (low) / 2 (intermediate) / 3 (high)  | 287     |
| `er`      | binary      | Estrogen receptor status: 0 negative / 1 positive        | 338     |
| `pr`      | binary      | Progesterone receptor status: 0 negative / 1 positive    | 356     |
| `status`  | binary      | 0 alive/censored, 1 died                                | 0       |
| `time`    | numeric     | Overall survival time (months)                           | 0       |

Outcome: overall survival (OS). Only 72 events (~6%), so this is a
low-event, low-dimensional survival setting — a deliberate stress-test for
comparing classical and machine-learning survival models.

Note: this is a public research dataset used for an educational project. If
you deploy or publish, verify its license terms on Kaggle.
