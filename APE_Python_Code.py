"""
Created on Fri Nov 16 12:51:19 2018

@author: Nguyen Thi Ngoc Mai
Student Number: 1746303
"""
#------------- Import necessary package for programming -----------------#
                                                                         #
from pyteomics import mgf, pepxml, mass                                  #
import numpy as np                                                       #
import matplotlib.pyplot as plt                                          #
import pandas as pd                                                      #
import os                                                                #
import function                                                          #
import pylab                                                             #
import math                                                              #
from statistics import mean , median                                     #
#------------------------------------------------------------------------#

##############################################################################
# ------------- STEP 1: Reading file to get mass and its charge ------------ #
##############################################################################

path = 'E:\\Study\\Uhasselt - Material\\2nd year\\Protein Analysis\\Project\\3800\\'
ms1 = dict()

for root, dirs, files in os.walk("3800"):  
    for filename in files:
        df = pd.read_csv(path+filename,sep = ' ', nrows = 1,skip_blank_lines=True, names = ['m/z','z'])
        ms1[df['m/z'][0]] = df['z'][0]


##############################################################################
# --------- STEP 2: Reading file to get mass and intensity from MS1 -------- #
##############################################################################

name = list(ms1.keys())
dataset = dict()

for root, dirs, files in os.walk("3800"):  
    for i in range(len(files)):
        #print(len(files))
        df = pd.read_csv(path+files[i], sep = ' ', skiprows = 1,
                         skip_blank_lines=True, names = ['m/z','intensity'])
        dataset[name[i]] = df

# Normalizing intensity of MS2
for i in name:
    df = dataset[i]
    df['i.n'] = [df['intensity'][j]/sum(df['intensity']) for j in range(len(df))]
    #df['r_intens'] = [df['intensity'][j]/max(df['intensity']) for j in range(len(df))]
    #df['n_intens'] = [df['intensity'][j]/mean(df['intensity']) for j in range(len(df))]
    #df['m_intens'] = [df['intensity'][j]/median(df['intensity']) for j in range(len(df))]
    #df['h'] = [df['intensity'][j]/df['m/z'][j] for j in range(len(df))]
    #df['norm_i'] = [df['h'][j]/sum(df['h']) for j in range(len(df))]
    # -- normalize intensity for observed data
     
##############################################################################
# -------- STEP 3: Creating database for protein, decoy and peptides --------#
##############################################################################
        
# ---------------- CREATE PROTEIN DATABASE FROM FASTA FILE ------------------#

path = 'E:\\Study\\Uhasselt - Material\\2nd year\\Protein Analysis\\Project\\'
p_protein = open(path+'studentP.fasta','r')

protein = p_protein.read().splitlines()

protein_lst = function.pro(protein)

protein_db = function.get_db(protein_lst)

# ----------------- CREATE DECOY DATABASE FROM FASTA FILE -------------------#

d_protein = open(path+'studentD.fasta','r')

decoy = d_protein.read().splitlines()

decoy_lst = function.pro(decoy)

decoy_db = function.get_decoydb(decoy_lst)


# ----------------- CREATE PEPTIDE DATABASE FROM FASTA FILE -----------------#

p_peptide = open(path+'studentP_peptides.fasta','r')
peptide = p_peptide.read().splitlines()

peptide_mass_db = function.get_mass_db(peptide)

d_peptide = open(path+'studentD_peptides.fasta','r')
decoy_p = d_peptide.read().splitlines()

decoy_mass_db = function.get_mass_db(decoy_p)


# --------------- CREATE DATABASE FOR PROTEIN AND ITS PEPTIDES----------------#

# 3rd party package to digest whole protein sequence to peptides
from pyteomics import parser

# merge protein and decoy database:
whole_db = dict()
for i in list(protein_db.keys()):
    whole_db[i] = protein_db[i]
for i in list(decoy_db.keys()):
    whole_db[i] = decoy_db[i]
    
protein_peptides = dict()
keys = list(whole_db.keys())

for i in keys:
    protein_peptides[i] = list(parser.cleave(whole_db[i], parser.expasy_rules['trypsin'],0))
        
    
##############################################################################
# -------- STEP 4: Getting potential candidates based on MS1's mass -------- #
##############################################################################

def candidate_ms1(precursor_mass, tolerance):
    '''
    This function to select potential candidate from current database
    of peptides digested from protein and decoy db
    - precursor_mass is from MS1 as MH+
    - tolerance value here is measured by ppm
    - proton: 1.007282 dalton (will be deducted from precursor_mass)
    '''
    protein = []
    decoy = []
    mass = precursor_mass - 1.00728 
    peptide_mass_lst = list(peptide_mass_db.keys())
    decoy_mass_lst = list(decoy_mass_db.keys())
    for i in peptide_mass_lst:      
        if abs(i - mass)/i*10**6 <= tolerance: #<= i + tolerance:
            protein.append(peptide_mass_db[i])
    for j in decoy_mass_lst:
        if abs(j - mass)/j*10**6 <= tolerance:
            decoy.append(decoy_mass_db[j])
    return protein, decoy

def matchPeptide2Protein(peptide):
    key = list(whole_db.keys())
    l = []
    for i in key:
        if peptide in protein_peptides[i]:
           l.append(i)
    return l

candidate=dict() # list of peptide sequences based on precursor mass
tolerance = 100 
for i in name:
    candidate[i]=candidate_ms1(i,tolerance)[0] + candidate_ms1(i,tolerance)[1]
    
ms1_candidate = {} # this is just to create dataframe for the output

for i in name:
    ms1_candidate[i] = pd.DataFrame([[j, matchPeptide2Protein(j)] for j in candidate[i]], columns = ['peptide', 'access.no'])
    #ms1_df[i] = pd.DataFrame(ms1_candidate[i],columns=['peptide', 'access.no'])

##############################################################################
# -------- STEP 5: Getting b,y ions from given sequence and its m/z -------- #
##############################################################################

def bIon_db(sequence,charge):
    '''
    Creat database for mass of b ions to compare with observed data
    '''
    b=function.bIon(sequence)
    b_db = dict()
    for j in range(1, charge):
        for i in b:
            b_db[float(mass.fast_mass(i,ion_type='b',charge=j))] = i# for j in range(1,charge)]
    return b_db
    
def yIon_db(sequence,charge):
    '''
    Creat database for mass of y ions to compare with observed data
    '''
    y=function.yIon(sequence)
    y_db = dict()
    for j in range(1, charge):
        for i in y:
            y_db[float(mass.fast_mass(i,ion_type='y',charge=j))] = i# for j in range(1,charge)]
    return y_db

def getIons(sequence,charge):
    '''
    This function return theoretical mass to charge of b and y ions
    from MS2.
    Based on this, comparing with observed data to list down
    potential ions to predict peptide sequence
    '''
    outcome = []
    bions = function.bIon(sequence)
    yions = function.yIon(sequence)
    for i in bions:
        outcome.append(i)
        for j in range(1,charge):
            outcome.append(float(mass.fast_mass(i,ion_type='b',charge=j)))
    for i in yions:
        outcome.append(i)
        for j in range(1,charge):
            outcome.append(float(mass.fast_mass(i,ion_type='y',charge=j)))
    return outcome

def getCIDFragmentIons(sequence,charge):
    l = getIons(sequence, charge)
    d = dict()
    out = []
    for i in range(int(len(l)/charge)):
        out.append(l[i*charge:i*charge+charge])
    for i in ('b', 'y'):
        d['b'] = pd.DataFrame(out[:int(len(out)/2)])
        d['y'] = pd.DataFrame(out[int(len(out)/2):])
    #return d
    for i in list(d.values()):
        print(i)
    
##### Example for function getCIDFragmentIons:
getCIDFragmentIons('AMLKPWTS',3)

# Results are matrices for b and y ions and their m/z with z up to 
# precursor charge - 1
    
#         0           1           2
#0        A   72.044386   36.525831
#1       AM  203.084876  102.046076
#2      AML  316.168936  158.588106
#3     AMLK  444.263896  222.635586
#4    AMLKP  541.316656  271.161966
#5   AMLKPW  727.395966  364.201621
#6  AMLKPWT  828.443646  414.725461
#         0           1           2
#0  STWPKLM  862.449131  431.728204
#1   STWPKL  731.408641  366.207959
#2    STWPK  618.324581  309.665929
#3     STWP  490.229621  245.618449
#4      STW  393.176861  197.092069
#5       ST  207.097551  104.052414
#6        S  106.049871   53.528574

def candidate_bIon(sequence, charge, ms1_mass, tolerance):
    '''
    This function to select potential b or y ions from potential candidate
    from MS1.
    It will return the list of b and y ions based on the mass in observed
    data within the tolerance
    Tolerance value is measured by da: this project applies 0.1da. User can 
    used by their own reference
    '''
    db = bIon_db(sequence,charge)
    ion = []
    #H = 1.00728
    for i in dataset[ms1_mass]['m/z']:
        for j in list(db.keys()):
            if abs(i-j)*charge <= tolerance:
                ion.append(db[j])
    return ion

def candidate_yIon(sequence, charge, ms1_mass, tolerance):
    '''
    This function to select potential b or y ions from potential candidate
    from MS1.
    It will return the list of b and y ions based on the mass in observed
    data within the tolerance
    Tolerance value is measured by da: this project applies 0.1da. User can 
    used by their own reference
    '''
    db = yIon_db(sequence,charge)
    ion = []
    for i in dataset[ms1_mass]['m/z']:
        for j in list(db.keys()):
            if abs(i-j)*charge <= tolerance:
                ion.append(db[j])
    return ion


##############################################################################
# ------- STEP 7: Scoring potential peptides from MS1 by count peak  ------- #
# --- by mass matching of ions and divided by number of theoretical ions --- #
##############################################################################

# Method 1: Counting number of peak matching
def score1(lst, precursor_mass):
    d = dict()
    for i in lst:  
        d[i] = len(candidate_bIon(i,ms1[precursor_mass],precursor_mass,0.1))+len(candidate_yIon(i,ms1[precursor_mass],precursor_mass,0.1))
    return d

# Method 2: Taking intensity and number of peak matching
# -- Getting the interval and average intensity within interval for observed data

def intensity_bIon(sequence, charge, ms1_mass, tolerance):
    '''
    This function to select potential b or y ions from potential candidate
    from MS1.
    It will return the list of b and y ions based on the mass in observed
    data within the tolerance
    Tolerance value is measured by da: this project applies 0.1da. User can 
    used by their own reference
    '''
    db = bIon_db(sequence,charge)
    ion = []
    s = 0
    df = dataset[ms1_mass]
    for i in range(len(df['m/z'])):
        for j in list(db.keys()):
            if abs(df['m/z'][i]-j)*charge <= tolerance:
                ion.append([db[j],df['m/z'][i]])
                s += df['i.n'][i]
    return s

def intensity_yIon(sequence, charge, ms1_mass, tolerance):
    '''
    This function to select potential b or y ions from potential candidate
    from MS1.
    It will return the list of b and y ions based on the mass in observed
    data within the tolerance
    Tolerance value is measured by da: this project applies 0.1da. User can 
    used by their own reference
    '''
    db = yIon_db(sequence,charge)
    ion = []
    s = 0
    df = dataset[ms1_mass]
    for i in range(len(df['m/z'])):
        for j in list(db.keys()):#Ion_db:
            if abs(df['m/z'][i]-j)*charge <= tolerance:
                ion.append([db[j],df['m/z'][i]])
                s += df['i.n'][i]
    return s

def score2(lst, precursor_mass):
    d = dict()
    for i in lst:  
        C = intensity_bIon(i,ms1[precursor_mass],precursor_mass,0.1)+intensity_yIon(i,ms1[precursor_mass],precursor_mass,0.1)
        f = math.exp(len(set(candidate_bIon(i,ms1[precursor_mass],precursor_mass,0.1)))+len(set(candidate_yIon(i,ms1[precursor_mass],precursor_mass,0.1))))
        d[i] = math.log10(C*f)
    return d
        
# Getting score for potential peptide candidates by each precursor mass
    
for i in name: # name is list of precursor mass
    print(pd.DataFrame(list(score2(candidate[i],i).values()),list(score2(candidate[i],i).keys()), columns=['Score for mass %s'%i]))

##############################################################################
# ------- STEP 8: Function for assigning confidence on the findings  ------- #
##############################################################################
def global_fdr(d):
    threshold = 0
    decoy = 0
    results = []
    for i in list(d.keys()):
        if d[i] > threshold:
            results.append(i)
    for j in results:
        if j in list(decoy_mass_db.values()):
            decoy += 1
    return round(decoy/len(results),4)

fdr=dict()
for i in name:
    fdr[i] = global_fdr(score2(candidate[i],i))
    
print(pd.DataFrame(list(fdr.values()),list(fdr.keys()), columns = ['FDR']))

# ------- List of potential peptide candidates in target and decoy db ------- #

target = dict()
decoy = dict()

for i in name:
    target[i] = candidate_ms1(i,100)[0]
    decoy[i] = candidate_ms1(i,100)[1]

