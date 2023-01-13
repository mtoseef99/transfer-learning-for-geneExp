AITL analysis
-------------

Download the processed data from this
[link](https://drive.google.com/drive/folders/1r4iw-qp6gep5XzYlQzaxI22cGgqI0MRG?usp=sharing),
provided by the authors, and create a [conda
environment](https://github.com/hosseinshn/AITL#conda-environment).

We have performed experiments for all four drugs (Bortezomib, Cisplatin,
Docetaxel, and Paclitaxel) with different k-folds cross-validation and
hyperparameter values given in \<these files\>.

 

These are the results for 24 experiments with baseline settings and
hyperparameter tuned models. The author performed all experiments with 3-fold
cross-validation, while we evaluated the model performance with also 5-fold and
10-fold cross validation.

| **Drug**   | **Model**            | **Cross-Validation** | **Epochs** | **Avg AUROC** | **Avg APR** |
|------------|----------------------|----------------------|------------|---------------|-------------|
| Bortezomib | Baseline             | 3-fold               | 15         | 0.732281153   | 0.742361278 |
| Bortezomib | Baseline             | 5-fold               | 15         | 0.7033083     | 0.727293207 |
| Bortezomib | Baseline             | 10-fold              | 15         | 0.737689394   | 0.737689394 |
| Cisplatin  | Baseline             | 3-fold               | 15         | 0.604796707   | 0.857139326 |
| Cisplatin  | Baseline             | 5-fold               | 15         | 0.622046784   | 0.85488424  |
| Cisplatin  | Baseline             | 10-fold              | 15         | 0.607037037   | 0.873774411 |
| Docetaxel  | Baseline             | 3-fold               | 15         | 0.492849743   | 0.62124888  |
| Docetaxel  | Baseline             | 5-fold               | 15         | 0.531868132   | 0.673401025 |
| Docetaxel  | Baseline             | 10-fold              | 15         | 0.441666667   | 0.647979025 |
| Paclitaxel | Baseline             | 3-fold               | 15         | 0.530762929   | 0.612060264 |
| Paclitaxel | Baseline             | 5-fold               | 15         | 0.515678084   | 0.58879782  |
| Paclitaxel | Baseline             | 10-fold              | 15         | 0.496616162   | 0.593677809 |
| Bortezomib | Hyperparameter tuned | 3-fold               | 15         | 0.720996487   | 0.732728673 |
| Bortezomib | Hyperparameter tuned | 5-fold               | 15         | 0.696646904   | 0.715587307 |
| Bortezomib | Hyperparameter tuned | 10-fold              | 15         | 0.732760295   | 0.751416034 |
| Cisplatin  | Hyperparameter tuned | 3-fold               | 15         | 0.55687164    | 0.823419209 |
| Cisplatin  | Hyperparameter tuned | 5-fold               | 15         | 0.604736842   | 0.844228739 |
| Cisplatin  | Hyperparameter tuned | 10-fold              | 15         | 0.578148148   | 0.852510261 |
| Docetaxel  | Hyperparameter tuned | 3-fold               | 15         | 0.498899249   | 0.647605641 |
| Docetaxel  | Hyperparameter tuned | 5-fold               | 15         | 0.522802198   | 0.666954559 |
| Docetaxel  | Hyperparameter tuned | 10-fold              | 15         | 0.448214286   | 0.654636072 |
| Paclitaxel | Hyperparameter tuned | 3-fold               | 15         | 0.550834613   | 0.608495571 |
| Paclitaxel | Hyperparameter tuned | 5-fold               | 15         | 0.546505152   | 0.62051767  |
| Paclitaxel | Hyperparameter tuned | 10-fold              | 15         | 0.546444444   | 0.618279505 |
