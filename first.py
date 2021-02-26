# -*- coding: utf-8 -*-
"""
Biomarkers

Created on Thu Jul 30 10:32:23 2020

@author: em322
"""

import csv
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from scipy.stats import pearsonr
import itertools as iter

def loadData(fName):
    '''
    Load data from csv file. Names in the first row,
    ids in the first column,
    data in other rows and columns
    
    return three arays:
        names
        ids
        data
    '''
    # Load data
    with open(fName, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        raw = []
        for row in reader:
            raw.append(row)
    # Organise arrays
    nCol = len(raw[0]) - 1
    names = raw[0][1:]
    raw = raw[1:]
    nRow = len(raw)
    ids = [None] * nRow
    data = np.zeros((nRow, nCol))
    for i in range(nRow):
        tmp = raw[i]
        ids[i] = tmp[0]
        for j in range(nCol):
            try:
                data[i, j] = float(tmp[j+1])
            except ValueError:
                data[i, j] = np.nan
    return names, ids, data

def testOneAttr(contr, ill, name, fout):
    '''
    Analised two samples contr and ill and write result to file fout.
    All NaNs are removed from ill an contr.
    if ill and contr contains less than 5 unique values then attribute is 
    considered as categorical. Otherwise attribute is tested as numerical.
    '''
    # Remove NaNs
    ill = ill[np.logical_not(np.isnan(ill))]
    contr = contr[np.logical_not(np.isnan(contr))]
    uniq = np.unique(np.concatenate((ill, contr)))
    # What do we have
    elems = uniq.shape[0]
    if elems == 0:
        print(name, 'is empty attribute\n', file=fout)
    elif elems == 1:
        print(name, 'is constant attribute\n', file=fout)
    elif elems < 5:
        # Convert data to int
        uniq = uniq.astype(np.int32)
        minVal = min(uniq)
        uniq = uniq - minVal
        ill = ill.astype(np.int32) - minVal
        contr = contr.astype(np.int32) - minVal
        maxValue = max(uniq) + 1
        cont = np.zeros((2, uniq.shape[0]))
        # Statistics for ill
        cnt = np.bincount(ill, minlength=maxValue)
        cnt = cnt[uniq]
        cont[1, :] = cnt
        percent = (cnt * (100 / sum(cnt))).astype(np.int32)
        cnt1 = cnt
        percent1 = percent
        # Statistics for contr
        cnt = np.bincount(contr, minlength=maxValue)
        cnt = cnt[uniq]
        cont[0, :] = cnt
        percent = (cnt * (100 / sum(cnt))).astype(np.int32)
        # Chi squared test
        st, p, dof, cnt2 = chi2_contingency(cont)
        st, pp, dof, cnt2 = chi2_contingency(cont, correction=False)
        print(name, '\tis categorical attribute\thealthy\t',\
              cnt, '\t', percent, '\till\t', cnt1, '\t', percent1, \
              '\tstat\t', pp, '\t', p, '\t', dof, '\n', file=fout)
        
    else:
        # Statistics for ill
        s, p = ttest_ind(ill, contr, equal_var = False)
        s, pp = mannwhitneyu(ill, contr, use_continuity=False, alternative='two-sided')
        print(name, '\tis numerical attribute\t',\
              'contr\t', contr.mean(), '\t', contr.std(),\
              '\till\t', ill.mean(), '\t', ill.std(),\
              '\tstat\t', p, '\t', pp, '\n', file=fout)

def pairedTestOneAttr(dat, name, fout):
    '''
    Analised one sample dat and write result to file fout.
    All NaNs are removed from dat.
    Considered only attributes with at least 5 different unique values. 
    All attributes are considered as numerical.
    '''
    # Remove NaNs
    dat = dat[np.logical_not(np.isnan(dat))]
    uniq = np.unique(dat)
    # What do we have
    elems = uniq.shape[0]
    if elems > 4:
        s, p = ttest_rel(dat, dat - dat)
        ss, pp = wilcoxon(dat)        
        print(name, '\t', dat.mean(), '\t', dat.std(), \
              '\t', p, '\t', pp, '\n', file=fout)

def formIndex(what, where):
    '''
    Form and return list of indices of elements of list what in the list where
    '''
    n = len(what)
    res = [None] * n
    for k in range(n):
        res[k] = where.index(what[k])
    return res

def removeNaN(rows, cols):
    ind = np.logical_not(np.isnan(rows))
    rows = rows[ind]
    cols = cols[ind]
    ind = np.logical_not(np.isnan(cols))
    rows = rows[ind]
    cols = cols[ind]
    return rows, cols
    

def calcCorrel(rows, cols):
    '''
    get variables from matrices rows and cols and calculate correlation matrix 
    with corr(rows[:, i],cols[:, j]) in element [i,j] and matrix of p-values
    '''
    nR = rows.shape[1]
    nC = cols.shape[1]
    res = np.zeros((nR, nC))
    pval = np.zeros((nR, nC))
    for r in range(nR):
        for c in range(nC):
            res[r, c], pval[r, c] = pearsonr(*removeNaN(rows[:, r], cols[:, c]))
    return(res, pval)
    
def writeTable(rows, cols, data, fout):
    '''
    write table with names of columns and rows in cols and rows and content in data.
    '''
    nC = len(cols)
    nR = len(rows)
    for c in range(nC):
        fout.write('\t{:s}'.format(cols[c]))
    for r in range(nR):
        fout.write('\n' + rows[r])
        for c in range(nC):
            fout.write('\t{:.2f}'.format(data[r, c]))
    fout.write('\n\n')
    
    
#
# General constants
#
# List of CRM outcomes
crmOuts = ["LV EDVi (mL?/m2)", "LV mass indexed to BSA (g?/m2)",\
           "LV mass?/volume (g?/ml)", "LV GLS (%)", "LV long PEDSR (1?/s)",\
           "LV circ PEDSR (1?/s)", "LV GCS (%)",\
           "Maxium LA volume indexed to BSA (ml?/m2)", "LA EF (%)",\
           "Average E?/e'", "LGE present (y?/n)", "Global MPR",\
           "Mean aortic distensibility (mmHg-1x10-3)"]
# List of biomarkers
biomarkers = ["Adiponectin Log10", "Angiopoietin 2 Log10", "BNP Log10",\
              "CRP BL Log10", "CHI3L Log10", "Cystatin C Log 10",\
              "Endoglin Log 10", "Endostatin Log 10", "ET1 Log 10",\
              "FABP3 Log 10", "FABP4 Log 10", "FGF21 Log 10", "FGF23 Log 10",\
              "Fas Log 10", "GDF15 Log 10", "Galectin 3 Log 10",\
              "ICAM1 Log 10", "IL1 beta Log 10", "IL10 Log 10", "IL6 Log 10",\
              "IL8 Log 10", "NGAL Log 10", "MMP12 Log 10", "MMP2 Log 10",\
              "MMP3 Log 10", "MMP7 Log 10", "MMP8 Log 10", "MMP9 Log 10",\
              "MPO Log 10", "NTproANP Log 10", "NTproBNP Log 10",\
              "OPG Log 10", "OPN Log 10", "p Selectin Log 10", "PAI1 Log 10",\
              "Pentraxin 3 Log 10", "Renin Log 10", "ST2 Log 10",\
              "Syndecan 1 Log 10", "Syndecan 4 Log 10", "KIM1 Log 10",\
              "TIMP1 Log 10", "TIMP4 Log 10", "TNFR1 Log 10", "TNFR2 Log 10",\
              "TNF alpha Log 10", "Tenascin C Log 10", "Troponin T Log 10",\
              "VEGFR1 Log 10", "VEGFa Log 10", "Leptin (pg?/L)"]

# Calculate individual statistics for attributes
attrStat = True
# Calculate correlation between variables
corrCalc = True

# To select one of the 5 "chapters"
# 1 is descriptive and correlation statistics for file Original.csv. Results will be written to file tests.txt
# 2 is descriptive statistics for subgroups in file Original.csv. Results will be written to file testsSubgr.txt
# 3 is descriptive statistics for subgroups in file After12weeks.csv. Results will be written to file testsSubgr12w.txt
# 4 is descriptive statistics for changes in subgroups in file Changes.csv. Results will be written to file testsSubgr12w.txt
# 5 is descriptive statistics for significance of changes in subgroups in file Changes.csv. Results will be written to file testsSubgrPair.txt

chapter = 5

if chapter == 1:
    #
    # Load data for work
    #
    names, ids, data = loadData("Original.csv")
    # remove columns with manu missing values in healthy group
    #ind2remove = [7, 15, 16, 17, 18, 19, 20, 75, 76, 89, 148, 149, 150]
    #names = np.delete(names, ind2remove).tolist()
    #data = np.delete(data, ind2remove, axis=1)
    
    # Separate two groups of patients
    nCol = len(names)
    contrId = ids[:36]
    contrData = data[:36,]
    illId = ids[36:]
    illData = data[36:,]
    
    if attrStat:
        # Prepare file for results
        with open("tests.txt", "w") as fout:
            # Get one attribute and send it to test
            for i in range(nCol):
                testOneAttr(contrData[:, i], illData[:, i], names[i], fout)
    
    if corrCalc:
        # Form vector of indices for crmOuts
        crmInd = formIndex(crmOuts, names)
        bioInd = formIndex(biomarkers, names)
        # Calculate correlation matrix for ill and healthy
        corrIll, pvalIll = calcCorrel(illData[:, bioInd], illData[:, crmInd])
        corrContr, pvalContr = calcCorrel(contrData[:, bioInd], contrData[:, crmInd])
        
        with open("corr.txt", "w") as fout:
            fout.write('ill\n')
            writeTable(np.asarray(names)[bioInd], np.asarray(names)[crmInd], corrIll, fout)
            fout.write('ill p-values\n')
            writeTable(np.asarray(names)[bioInd], np.asarray(names)[crmInd], pvalIll, fout)
            fout.write('Corr\n')
            writeTable(np.asarray(names)[bioInd], np.asarray(names)[crmInd], corrContr, fout)
            fout.write('Corr p-values\n')
            writeTable(np.asarray(names)[bioInd], np.asarray(names)[crmInd], pvalContr, fout)
            fout.write('Differences\n')
            writeTable(np.asarray(names)[bioInd], np.asarray(names)[crmInd], corrContr - corrIll, fout)
        
if chapter == 2:
    #
    # Load data for work
    #
    names, ids, data = loadData("Original.csv")
    # Now compare three groups of ill participants

    # Firstly remove attributes with many missing values. I decided to not 
    # remove because of each column has at least 11 values and some coparison
    # can be done.   
    
    # Split data into three groups
    # Get Subgroup column
    ind = formIndex(['Subgroup'], names)
    subgr = (data[:, ind]).flatten()
    # Extract subgroups 2, 3, 4
    ind = subgr == 2
    data2 = data[ind, :]
    id2 = list(iter.compress(ids,ind))
    ind = subgr == 3
    data3 = data[ind, :]
    id3 = list(iter.compress(ids,ind))
    ind = subgr == 4
    data4 = data[ind, :]
    id4 = list(iter.compress(ids,ind))
    
    # Prepare file for results
    if attrStat:
        nCol = len(names)
        with open("testsSubgr.txt", "w") as fout:
            # Get one attribute and send it to test
            print('Comparison of 2 and 3', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data3[:, i], names[i], fout)
            print('\n\n\nComparison of 2 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data4[:, i], names[i], fout)
            print('\n\n\nComparison of 3 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data3[:, i], data4[:, i], names[i], fout)
    
if chapter == 3:
    #
    # Load data for work
    #
    names, ids, data = loadData("After12weeks.csv")
    # Now compare three groups of ill participants

    # Firstly remove attributes with many missing values. I decided to not 
    # remove because of each column has at least 11 values and some coparison
    # can be done.   
    
    # Split data into three groups
    # Get Subgroup column
    ind = formIndex(['Treatment Group'], names)
    subgr = (data[:, ind]).flatten()
    # Extract subgroups 2, 3, 4
    ind = subgr == 2
    data2 = data[ind, :]
    id2 = list(iter.compress(ids,ind))
    ind = subgr == 3
    data3 = data[ind, :]
    id3 = list(iter.compress(ids,ind))
    ind = subgr == 4
    data4 = data[ind, :]
    id4 = list(iter.compress(ids,ind))
    
    # Prepare file for results
    if attrStat:
        nCol = len(names)
        with open("testsSubgr12w.txt", "w") as fout:
            # Get one attribute and send it to test
            print('Comparison of 2 and 3', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data3[:, i], names[i], fout)
            print('\n\n\nComparison of 2 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data4[:, i], names[i], fout)
            print('\n\n\nComparison of 3 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data3[:, i], data4[:, i], names[i], fout)
    
if chapter == 4:
    #
    # Load data for work
    #
    names, ids, data = loadData("Changes.csv")
    # Now compare three groups of ill participants

    # Firstly remove attributes with many missing values. I decided to not 
    # remove because of each column has at least 11 values and some coparison
    # can be done.   
    
    # Split data into three groups
    # Get Subgroup column
    ind = formIndex(['Treatment Group'], names)
    subgr = (data[:, ind]).flatten()
    # Extract subgroups 2, 3, 4
    ind = subgr == 2
    data2 = data[ind, :]
    id2 = list(iter.compress(ids,ind))
    ind = subgr == 3
    data3 = data[ind, :]
    id3 = list(iter.compress(ids,ind))
    ind = subgr == 4
    data4 = data[ind, :]
    id4 = list(iter.compress(ids,ind))
    
    # Prepare file for results
    if attrStat:
        nCol = len(names)
        with open("testsSubgrCh.txt", "w") as fout:
            # Get one attribute and send it to test
            print('Comparison of 2 and 3', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data3[:, i], names[i], fout)
            print('\n\n\nComparison of 2 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data2[:, i], data4[:, i], names[i], fout)
            print('\n\n\nComparison of 3 and 4', file=fout)
            for i in range(nCol):
                testOneAttr(data3[:, i], data4[:, i], names[i], fout)
    
if chapter == 5:
    #
    # Load data for work
    #
    names, ids, data = loadData("Changes.csv")
    # Now compare three groups of ill participants

    # Firstly remove attributes with many missing values. I decided to not 
    # remove because of each column has at least 11 values and some coparison
    # can be done.   
    
    # Split data into three groups
    # Get Subgroup column
    ind = formIndex(['Treatment Group'], names)
    subgr = (data[:, ind]).flatten()
    # Extract subgroups 2, 3, 4
    ind = subgr == 2
    data2 = data[ind, :]
    id2 = list(iter.compress(ids,ind))
    ind = subgr == 3
    data3 = data[ind, :]
    id3 = list(iter.compress(ids,ind))
    ind = subgr == 4
    data4 = data[ind, :]
    id4 = list(iter.compress(ids,ind))
    
    # Prepare file for results
    if attrStat:
        nCol = len(names)
        with open("testsSubgrPair.txt", "w") as fout:
            # Get one attribute and send it to test
            print('Group 2', file=fout)
            for i in range(nCol):
                pairedTestOneAttr(data2[:, i], names[i], fout)
            print('Group 3', file=fout)
            for i in range(nCol):
                pairedTestOneAttr(data3[:, i], names[i], fout)
            print('Group 4', file=fout)
            for i in range(nCol):
                pairedTestOneAttr(data4[:, i], names[i], fout)
