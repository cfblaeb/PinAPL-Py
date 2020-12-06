#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:22:56 2016

@author: philipp
"""
# Library sanity check
# =======================================================================
# Imports
import os
import pandas


def RunSanityCheck(config):
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    LibDir = config['LibDir']
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    libsep = '\t' if LibFormat == 'tsv' else ','
    DataDir = config['DataDir']
    WorkingDir = config['WorkingDir']
      
    # --------------------------------------------------------------------
    # Replace non-printable characters from library (...these cause problems in PlotCount.py)
    # --------------------------------------------------------------------
    LibCols = ['gene', 'ID', 'seq']
    LibFile = pandas.read_csv(f"{LibDir}/{LibFilename}", sep=libsep, skiprows=1, names=LibCols)
    GeneNames = list(LibFile['gene'])
    ID = list(LibFile['ID'])
    seq = list(LibFile['seq']) 
    GeneNames0 = []
    ID0 = []
    
    # --------------------------------------------------------------------
    # Define bad characters (library)
    # --------------------------------------------------------------------     
    BadCharacters = [' ','>','<',';',':',',','|','/','\\','(',')','[',']', '$','%','*','?','{','}','=','+','@']
        
    # --------------------------------------------------------------------
    # Check library
    # --------------------------------------------------------------------         
    BadLibCharFound = False
    for gene in GeneNames:
        for bad_char in BadCharacters:
            gene = gene.replace(bad_char,'_')
        GeneNames0.append(gene)
    for sgRNA in ID:
        for bad_char in BadCharacters:
            sgRNA = sgRNA.replace(bad_char,'_')       
        ID0.append(sgRNA)    
    if GeneNames != GeneNames0 or ID != ID0:
        BadLibCharFound = True
        LibFile0 = pandas.DataFrame(data = {'gene': [gene for gene in GeneNames0], 'ID': [sgRNA for sgRNA in ID0], 'seq': [s for s in seq]}, columns = ['gene','ID','seq'])
        LibFile0.to_csv(f"{LibDir}/LibFilename", sep = libsep, index = False)
        print("WARNING: Special characters in library file have been replaced by '_' ")

    # --------------------------------------------------------------------
    # Load Data Sheet
    # --------------------------------------------------------------------
    DataSheet = pandas.read_excel(f'{WorkingDir}/DataSheet.xlsx')
    Filenames = list(DataSheet['FILENAME'])
    TreatmentList = list(DataSheet['TREATMENT'])

    # --------------------------------------------------------------------
    # Define bad characters (filenames & samples)
    # --------------------------------------------------------------------         
    BadCharacters = [' ','>','<',';',':',',','|','/','\\','(',')','[',']', '$','%','*','?','{','}','=','+','@']
    
    # --------------------------------------------------------------------
    # Replace non-printable characters from filenames 
    # --------------------------------------------------------------------
    BadFileCharFound = False        
    for j in range(len(Filenames)):
        Filename = Filenames[j]
        Filename0 = Filename
        for bad_char in BadCharacters:
            Filename0 = Filename0.replace(bad_char,'_')
        if Filename0 != Filename:
            BadFileCharFound = True
            os.system(f"mv {DataDir}/{Filename} {DataDir}/{Filename0}")
            DataSheet['FILENAME'][j] = Filename0  
            print("WARNING: Special characters in filenames names replaced by '_'")

    # --------------------------------------------------------------------
    # Replace non-printable characters from sample names 
    # --------------------------------------------------------------------   
    TreatmentList0 = TreatmentList
    BadSampleCharFound = False 
    for bad_char in BadCharacters:
        TreatmentList0 = [str(treatment).replace(bad_char,'_') for treatment in TreatmentList0]
    if TreatmentList0 != TreatmentList:
        BadSampleCharFound = True
        DataSheet['TREATMENT'] = TreatmentList0        
        print("WARNING: Special characters in sample names replaced by '_'")
        
    # --------------------------------------------------------------------
    # Update Data Sheet
    # --------------------------------------------------------------------            
    if BadFileCharFound or BadSampleCharFound:
        DataSheet.to_excel(f'{WorkingDir}/DataSheet.xlsx',columns=['FILENAME','TREATMENT'])

    # --------------------------------------------------------------------
    # No special characters found
    # -------------------------------------------------------------------- 
    if not BadLibCharFound and not BadFileCharFound and not BadSampleCharFound:
        print('No special characters found.')
