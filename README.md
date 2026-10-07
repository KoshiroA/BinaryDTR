This repository provides a reproduction program for the simulations presented in:

"Regression-based doubly robust estimation of optimal dynamic treatment regimes for binary outcomes using collapsible effect measures"

by Koshiro Arai and Tomohiro Shinozaki

History:

⦁	20260124:　First Upload
⦁	20261007:　Second Upload

---
[0.Preparation]: Library and Estimation method

**"0.1_Libs.R"**

**"0.2_QLearn.R"**: Performing Q-learning

**"0.3_Gestimation.R"**: Performing G-estimation with logistic nuisance treatment-free model

**"0.6_Gest_OP.R"**: Performing G-estimation with logistic odds product model
  
  Before running **"0.6_Gest_OP.R"**, perform them in order:
  **"0.4_Qk_estimation.R"**; **"0.5_Solve_Gestimation_equation.R"**

---

[Simulation 1]

**"1.1_Datagen_logi.R"**

**"1.2_VsimFunc_logi.R"**

**"1.3_Run_Sim_logi.R"**　-> Table 2, Figure 2; Table S1, Figure S1, Figure S4

---
[Simulation 2]

**"2.1_Datagen_opt_general.R"**

**"2.2_VsimFunc_opt.R"**

**"2.3_Run_Sim_opt.R"**-> Table 3, Figure 3,4; Table S1, Figure S1, Table S2, Table S3, Figure S2-S3, S5-S6

[Additional Simulation in Section S2.4]

**"3.1_Datagen_DR.R"**

**"3.2_Gmod_DR.R"**

**"3.3_GSim_DR.R"**

**"3.3_GSim_DR.R"**-> Table S6

---
[Seed_summary]

Seed_(Sim Number)(Effect scenarios)_(SampleSize)

ex. Seed_1a_5000: Random seed after running the Simulation 1 of n = 5000 under alternative hypothesis

ex. Seed_2a_500: Random seed after running the Simulation of 2 n = 5000 under (a): ψk=(0.3,-0.45)
