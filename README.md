#  Regression-based doubly robust estimation of optimal dynamic treatment regime for binary outcomes using collapsible effect measures #

by Koshiro Arai and Tomohiro Shinozaki

⦁	20260124:　First Upload

---
[0.Preparation]
Library and Estimation method

0.1_Libs.R
0.2_QLearn.R
0.3_Gestimation.R
0.4_Gest_OP.R

To obtain doubly robust g-estimators, perform the code in order

--------Simulation in the main text and supplementary material--------

[Simulation 1]

1.1_Datagen_logi
1.2_VsimFunc_logi
1.3_Run_Sim_logi.R　-> Table 2, Figure 2; Table S1, Figure S1

[Simulation 2]

2.1_VsimFunc_opt
2.2_VsimFunc_opt
2.3_Run_Sim_opt.R-> Table 3, Figure 3,4; Table S1, Figure S1, Table S2, Table S3, Figure S2-S4

[Seed_summary]

Seed_(Sim Number)(Effect scenarios)_(SampleSize)

ex. Seed_1a_5000: Random seed after running the Simulation 1 of n = 5000 under alternative hypothesis
ex. Seed_2a_500: Random seed after running the Simulation of 2 n = 5000 under (a): ψk=(0.3,-0.45)
