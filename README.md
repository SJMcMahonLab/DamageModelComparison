# DamageModelComparison

TOPAS-nBio 1.0 (TOPAS 3.6.1) parameter files, MATLAB R2023a codes and Medras references that have been used to simulate Monte Carlo models of radiation-induced Double Strand Break (DSB) damage (Schuemann et al 2019a). This was completed in the Standard DNA Damage (SDD) data format to be used in combination with the Medras biological response model (McMahon and Prise 2021), to predict the yields of lethal chromosome aberrations. Damage models were designed with various levels of detail, to investigate the impact of DSB yield. Models are listed by reducing levels of simulation detail investigating the inclusion of chemical interactions, realistic nuclear geometry, single strand break damage and track structure. Each damage model underwent the same repair modelling to determine any impact on chromosome aberration yields.

*ChemistryModel.txt*

Simulates chemistry, full nuclear geometry, single strand break damage and track structure.
This TOPAS-nBio parameter file simulates the physical and chemical irradiation of a human fibroblast nucleus (Zhu et al 2020). Default direct and indirect SB damage, and scavenging parameters are specified but can be updated alongside the exposure conditions.

*PhysicsOnlyModel.txt*

Simulates full nuclear geometry, single strand break damage and track structure.
Same as the ChemistryModel.txt but with the simulation of chemical interactions removed and the direct damage parameter updated.

*Support Files*

Contain nuclear geometry and chemistry information to be read in by "ChemistryModel.txt" and "PhysicsOnlyModel.txt". Support files are sourced and available from TOPAS-nBio (Schuemann et al 2019b, [1]).

*SimpleDNAModel.m*

Simulates single strand break damage and track structure.
This MATLAB file reads in tuple scoring of physical interaction information from "ChemistryModel.txt" or "PhysicsOnlyModel.txt". 

*InitialDSBModel.m*

Simulates track structure.
This MATLAB file reads in tuple scoring of physical interaction information from "ChemistryModel.txt" or "PhysicsOnlyModel.txt". 

*MedrasDamageModel*

This is the most simple damage used and is available and detailed elsewhere (McMahon and Prise 2021, [2])

*chromModel.py*

Sourced from [2], this assigns chromosome positions to the DSB damage site in the Simple DNA, Initial DSB and Medras models for repair simulation.

*Convert_to_SDD*

Converts the Simple DNA and Initial DSB damage models DSB outputs to SDD format and assigns chromosome positions using "chromModel.py".

*MedrasRepairModel*

This model is available and detailed elsewhere (McMahon and Prise 2021, [2]). It inputs the SDD data files and simulates repair and response, giving misrepair rates and yields of aberrations.

# Contact
For any questions or problems regarding the codes, please contact sthompson@qub.ac.uk.

# References and links
[1] https://github.com/topas-nbio

[2] https://github.com/sjmcmahon/Medras-MC

McMahon S J and Prise K M 2021 A Mechanistic DNA Repair and Survival Model (Medras): Applications to Intrinsic Radiosensitivity, Relative Biological Effectiveness and Dose-Rate Front Oncol 11

Schuemann J, McNamara A L, Ramos-Méndez J, Perl J, Held K D, Paganetti H, Incerti S and Faddegon B 2019a TOPAS-nBio: An Extension to the TOPAS Simulation Toolkit for Cellular and Sub-cellular Radiobiology Radiat Res 191 125–38 

Schuemann J, Mcnamara A L, Warmenhoven J W, Henthorn N T, Kirkby K J, Merchant M J, Ingram S, Paganetti H, Held K D, Ramos-Mendez J, Faddegon B, Perl J, Goodhead D T, Plante I, Rabus H, Nettelbeck H, Friedland W, Kundrát P, Ottolenghi A, Baiocco G, Barbieri S, Dingfelder M, Incerti S, Villagrasa C, Bueno M, Bernal M A, Guatelli S, Sakata D, Brown J M C, Francis Z, Kyriakou I, Lampe N, Ballarini F, Carante M P, Davídková M, Štepán V, Jia X, Cucinotta F A, Schulte R, Stewart R D, Carlson D J, Galer S, Kuncic Z, Lacombe S, Milligan J, Cho S H, Sawakuchi G, Inaniwa T, Sato T, Li W, Solov’yov A V., Surdutovich E, Durante M, Prise K M and Mcmahon S J 2019b A New Standard DNA Damage (SDD) Data Format Radiat Res 191 76–92 

Zhu H, McNamara A L, McMahon S J, Ramos-Mendez J, Henthorn N T, Faddegon B, Held K D, Perl J, Li J, Paganetti H and Schuemann J 2020 Cellular Response to Proton Irradiation: A Simulation Study with TOPAS-nBio Radiat Res 194 9–21
