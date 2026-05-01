# emorl_paper_phaneuf-hadd

Data and Code

*to accompany*

Phaneuf-Hadd, C.V., Phelps, E.A., & Somerville, L.H. (in preparation). Changing Contributions of Emotion and Reward Information to Learning Through Adolescence.

## Developer Contact Information

Github profile: https://github.com/cphaneuf

Email: (current) cphaneuf@g.harvard.edu, (permanent) cphaneuf@umich.edu

## Contents

### data/ directory

Contains demographic, behavioral (learning task, test phase, scale), and computational modeling (model comparison, model recovery, parameter recovery, parameter fit) data.

### analyses/ directory

*utilities.R* defines variables and functions to be shared across scripts.

*demog.R* takes demographic and WASI matrix reasoning data inputs from data/ and data/behav/ and writes outputs to results/demog/.

*stat_mods.R* takes learning task, test phase, pre- and post-subjective rating, reinforcement report, and subjective value of money data inputs from data/behav/ and writes outputs to results/learn/, results/test/, and results/scales/.

*comp_mods.R* takes model comparison, model recovery, parameter recovery, and parameter fit secondary data (derived from primary learning task data) inputs from data/model_comp/, data/model_recov/, data/param_recov/, and data/param_fits/ and writes outputs to results/learn/\<subdirectories with the same names as in data>.

### results/ directory

Contains text and png file outputs from scripts in analyses/, sorted by analysis type.

### annotated_figs/ directory

Contains annotated_figs.pptx, which annotates several figures beyond the limits of R.

