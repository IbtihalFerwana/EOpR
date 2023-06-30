

EOpR
=======================================

Code for manuscript **Optimal Recovery for Causal Inference** https://arxiv.org/abs/2208.06729


Empirical Experiments Replication
========================

#### Abadie's Basque Country
```
python python/main.py --data data/basque_abadi_reshaped.csv --pretreatment_end 1971 --results results_basq --figures figures_basq
```
#### Abadie's California Proposition 99
```
python python/main.py --data data/prop99_abadie_reshaped.csv --pretreatment_end 18 --results results_cal --figures figs_cal
```

Citing
==========================
*BibTeX*:: 
```
  @misc{FerwanaV2022,
  doi = {10.48550/arXiv.2208.06729},
  url = {https://arxiv.org/abs/2208.06729},
  author = {Ferwana, Ibtihal and Varshney, Lav R.},
  keywords = {Methodology (stat.ME); Econometrics (econ.EM); Signal Processing (eess.SP)},
  title = {Optimal Recovery for Causal Inference},
  publisher = {arXiv},
  year = {2022},
  copyright = {Creative Commons Attribution 4.0 International}}
```
