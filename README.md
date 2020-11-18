License
=========

Copyright (C) 2020 Guihua Duan(duangh@mail.csu.edu.cn),Lingzhi Zhu(lz_zhu@csu.edu.cn)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.

Guihua Duan(duangh@mail.csu.edu.cn),Lingzhi Zhu(lz_zhu@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


Prediction of microbe-drug associations based on chemical structures and KATZ measure
=================
HMDAKATZ is one novel computational model, which utilizes KATZ to identify potential microbe-drug associations. The codes of HMDAKATZ are implemented in Matlab2015b.

1.DATA.

1.1 Drug_Number.xlsx store drug ids and names;

1.2 Microbe_Number.xlsx store microbe ids and names;

1.3 Microbe-Drug_Number.xlsx store the known microbe-drug associations;


2.HMDAKATZ_5CV.

2.1 knowndrugmicrobeinteraction.mat store the known microbe-drug associations;

2.2 drugsmiles.mat store the SMILES similarity of drugs;

2.3 main.m : 5-fold cross validation of main function;

2.4 HMDAKATZ_5CV.m: function of 5-fold cross validation of HMDAKATZ; 

2.5 positiontooverallauc.m: function of calculating the AUC value; 

3.HMDAKATZ_10CV.

3.1 knowndrugmicrobeinteraction.mat store the known microbe-drug associations;

3.2 drugsmiles.mat store the SMILES similarity of drugs;

3.3 main.m : 10-fold cross validation of main function;

3.4 HMDAKATZ_10CV.m: function of 10-fold cross validation of HMDAKATZ; 

3.5 positiontooverallauc.m: function of calculating the AUC value; 

4.HMDAKATZ_LOOCV.

4.1 knowndrugmicrobeinteraction.mat store the known microbe-drug associations;

4.2 drugsmiles.mat store the SMILES similarity of drugs;

4.3 main.m : Leave one cross validation of main function;

4.4 HMDAKATZ_LOOCV.m: function of leave one cross validation of HMDAKATZ; 

4.5 positiontooverallauc.m: function of calculating the AUC value; 

All files of Data and Code should be stored in the same folder to run HMDAKATZ.

Predict HMDA based on chemical structures and KATZ measure
