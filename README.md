# DamageModelComparison

TOPAS-nBio 1.0 (TOPAS 3.6.1) and Medras (verison 15/08/2023) parameter files and references (Schuemann et al 2019a, McMahon and Prise 2021), alongside additional MATLAB (R2023a) and Python analysis codes. These simulation toolkits and additional codes have been used to model radiation-induced Double Strand Break (DSB) damage in a nuclear volume. This was completed in the Standard DNA Damage (SDD) data format (Schuemann et al 2019b). DSB damage outputs were imported into the Medras biological response model to predict the yields of lethal chromosome aberrations. 

Damage models were designed with various levels of detail to investigate the impact of model design on the initial DSB damage and the formation of lethal aberrations. Models are listed by reducing levels of simulation detail investigating the inclusion of chemical interactions, realistic nuclear geometry, single strand break damage and track structure. 

*ChemistryModel.txt*

Simulates chemistry, full nuclear geometry, single strand break damage and track structure.
This TOPAS-nBio parameter file simulates the physical and chemical interactions during the irradiation of a human fibroblast nucleus (Zhu et al 2020). Default direct and indirect SB damage, and scavenging parameters are specified but can be updated alongside the exposure conditions.

*PhysicsOnlyModel.txt*

Simulates full nuclear geometry, single strand break damage and track structure.
Same as the ChemistryModel.txt but with the simulation of chemical interactions removed and the direct damage parameter updated.

*Support Files for Chemistry and Physics Only models*

The "ChemistryModel.txt" and "PhysicsOnlyModel.txt read in support files containing the nuclear geometry and chemistry. Support files are sourced and available from TOPAS-nBio (Schuemann et al 2019b, [1]).

*SimpleDNAModel.m*

Simulates single strand break damage and track structure.
This MATLAB file reads in tuple scoring of physical interaction information from "ChemistryModel.txt" or "PhysicsOnlyModel.txt". 

*InitialDSBModel.m*

Simulates track structure.
This MATLAB file reads in tuple scoring of physical interaction information from "ChemistryModel.txt" or "PhysicsOnlyModel.txt". 

*Medras Damage Model*

This is the most simple damage model used and is available and detailed elsewhere (McMahon and Prise 2021, [2])

*chromModel.py*

Sourced from [2], this assigns chromosome positions to the DSB damage site in the Simple DNA, Initial DSB and Medras models for repair simulation.

*Convert_to_SDD.py*

Converts the Simple DNA and Initial DSB damage models DSB outputs to SDD format and assigns chromosome positions using "chromModel.py".

*Medras Repair Model*

This model is available and detailed elsewhere (McMahon and Prise 2021, [2]). The SDD files from each damage model were input into Medras, simulating the repair and biological response to DSB damage, outputting misrepair rates and yields of aberrations.

# Contact
For any questions or problems regarding the codes, please contact sthompson@qub.ac.uk.

# References and links
[1] https://github.com/topas-nbio/TOPAS-nBio/tree/master/examples/geometry/nucleusModel/supportFiles

[2] https://github.com/sjmcmahon/Medras-MC

McMahon S J and Prise K M 2021 A Mechanistic DNA Repair and Survival Model (Medras): Applications to Intrinsic Radiosensitivity, Relative Biological Effectiveness and Dose-Rate Front Oncol 11

Schuemann J, McNamara A L, Ramos-Méndez J, Perl J, Held K D, Paganetti H, Incerti S and Faddegon B 2019a TOPAS-nBio: An Extension to the TOPAS Simulation Toolkit for Cellular and Sub-cellular Radiobiology Radiat Res 191 125–38 

Schuemann J, Mcnamara A L, Warmenhoven J W, Henthorn N T, Kirkby K J, Merchant M J, Ingram S, Paganetti H, Held K D, Ramos-Mendez J, Faddegon B, Perl J, Goodhead D T, Plante I, Rabus H, Nettelbeck H, Friedland W, Kundrát P, Ottolenghi A, Baiocco G, Barbieri S, Dingfelder M, Incerti S, Villagrasa C, Bueno M, Bernal M A, Guatelli S, Sakata D, Brown J M C, Francis Z, Kyriakou I, Lampe N, Ballarini F, Carante M P, Davídková M, Štepán V, Jia X, Cucinotta F A, Schulte R, Stewart R D, Carlson D J, Galer S, Kuncic Z, Lacombe S, Milligan J, Cho S H, Sawakuchi G, Inaniwa T, Sato T, Li W, Solov’yov A V., Surdutovich E, Durante M, Prise K M and Mcmahon S J 2019b A New Standard DNA Damage (SDD) Data Format Radiat Res 191 76–92 

Zhu H, McNamara A L, McMahon S J, Ramos-Mendez J, Henthorn N T, Faddegon B, Held K D, Perl J, Li J, Paganetti H and Schuemann J 2020 Cellular Response to Proton Irradiation: A Simulation Study with TOPAS-nBio Radiat Res 194 9–21
