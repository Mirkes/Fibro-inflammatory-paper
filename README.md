# Fibro-inflammatory-paper
This repository contains all data files and all code used for the Fibro-inflammatory-paper. List of authors and form of citation will be added after the paper publishing.

This repository contains several functions developed for this problem, several scripts for data processing. Scripts also used three more repositories:

Elastic maps https://github.com/Mirkes/ElMap used for drawing of points in the space of the first three principal components.

KNN based data imputation https://github.com/Mirkes/DataImputation for the data imputation to use in visualisation.

Separability dimension analysis https://github.com/auranic/FisherSeparabilityAnalysis/tree/master/MATLAB for dimensionality analysis.

## Python part of project - descriptive statistics
File first.py contains python code of descriptive statistics and correlation analysis. This file uses three csv files with data:

Original.csv contains data for the beginning.

After12weeks.csv contains data after 12 weeks.

Changes.csv contains changes between the beginning and state after 12 weeks.

## Matlab part of project
Below you can find description of Matlab functions and scripts.

### Functions:
1. oneDClass.m applied classification with one input attribute by searching the best threshold. Values that greater than or equal to threshold belongs to one class and Values that are less than threshold belongs to another class. This function calculate the best threshold and corresponding error. The used error measure is 1 - 0.5*(TP/Pos+TN/Neg), where 

&nbsp;&nbsp;&nbsp;TP is true positive or the number of correctly recognised classes of the first class, 
	
Pos is the number of classes of the first class,
	
	TN is true negative or the number of correctly recognised classes of the second class, 
	
	Neg is the number of classes of the second class. 
	
2. fisher.m This function calculate direction of fisher discriminant and then used oneDClass to find the best threshold and corresponding error. This function also return fisher discriminant direction.

4. saveFigures simple save figures to specified folder.

### Scripts:
a_Bio.m Form figures in the space of demographic variables and biomarkers variables. Used Bio.mat.

a_CMR.m Form figures in the space of demographic variables and CMR variables. Used CMR.mat.

a_Total.m Form figures in the space of demographic variables, CMR variables and biomerkers variables. Used Total.mat.

Changes.m Form figures for changes of points lokaction before and after 12 weeks in principal components. Used Changes.mat or ChangesFiltered.mat.

centreMoves.m Form figures in the plain of changes. Used ChangesFiltered.mat.

OneDand2D.m form figures for the best one attribute classifiers for each of subgroups; Search the best set of two attributes (the best Fisher discriminant separation), Calculate the full dimensional Fisher’s discriminant, applied dimensionality analysis, draw graph of PCs informativeness. Used Changes.mat or ChangesFiltered.mat.

### Data sets:
Bio.mat is set of demographic and biomarkers data. Fime contains matrix of data (BioData), list of groups 0 for healthy volunteers, 1 for standard care, 2 for exercise, and 3 for MRP (BioGroups), identifiers of participants (BioIDs), and list of attributes (Attrs). 

CMR.mat is set of demographic and CMR data. Fime contains matrix of data (CMRData), list of groups 0 for healthy volunteers, 1 for standard care, 2 for exercise, and 3 for MRP (CMRGroups), identifiers of participants (CMRIDs), and list of attributes (Attrs). 

Total.mat is set of demographic, CMR and biomarkers data. Fime contains matrix of data (TotalData), list of groups 0 for healthy volunteers, 1 for standard care, 2 for exercise, and 3 for MRP (TotalGroups), identifiers of participants (TotalIDs), and list of attributes (Attrs). 

Changes.mat contains data before start (Original) and after 12 weeks (After). Ids contains identifiers of patients and Columns contains names of attributes.

ChangesFiltered.mat contains data of healthy volunteers (Healthy), before start (Original) and after 12 weeks (After) with reduced set of attributes. Ids contains identifiers of patients and Columns contains names of attributes.
