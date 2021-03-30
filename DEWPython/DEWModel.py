#!/usr/bin/env python
# coding: utf-8

# # Package Imports 

# In[2]:


import pandas as pd
import numpy as np
import sys
import threading
import subprocess
import os
import json
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from DEWPython import DEWEquations
from matplotlib.lines import Line2D
get_ipython().run_line_magic('matplotlib', 'inline')
from collections import defaultdict
import pkg_resources
import os.path as op

# ### Defining a Global Variables (Location and Constants)
resource_package = 'DEWPython'
min_path = '/'.join(('resources', 'mineralDictionary.txt'))
mineralPath = pkg_resources.resource_filename(resource_package, min_path)

min_path2 = '/'.join(('resources', 'extendMineralDictionary.txt'))
mineralPath2 = pkg_resources.resource_filename(resource_package, min_path2)

gas_path = '/'.join(('resources', 'gasLst.txt'))
gasPath = pkg_resources.resource_filename(resource_package, gas_path)
aq_path = '/'.join(('resources', 'aqueousLst.txt'))
aqPath = pkg_resources.resource_filename(resource_package, aq_path)
die_path =  '/'.join(('resources', 'dielectric.csv'))
diePath = pkg_resources.resource_filename(resource_package, die_path)
inp_path ='/'.join(('resources', 'input.csv'))
inpPath = pkg_resources.resource_filename(resource_package, inp_path)
den_path ='/'.join(('resources', 'input.csv'))
denPath = pkg_resources.resource_filename(resource_package, den_path)
g_path = '/'.join(('resources', 'water_gibbs.csv'))
gPath = pkg_resources.resource_filename(resource_package, g_path)
sup_path = '/'.join(('resources', 'supcrt96.exe'))
sup_Path =  pkg_resources.resource_filename(resource_package, sup_path)

global Tr, bigQ, Chi, Pr, E_PrTr, bigR, Psi, Theta, Upsilon, Conversion, mineralDictionary

mineralDictionary = json.load(open(mineralPath))
mineralDictionary2 = json.load(open(mineralPath2))
'''A dictionary that stores all the minerals and allows them to be queried for use in the DEW model.'''

bigQ = 5.903E-07
'''Big Q is the 5.903E-07, and has units of bar^-1 '''
Chi = -3.090E-07
'''X is the constant -3.090E-07 and has units of K^-2'''
T_r = 298.15
'''The standard state temperature 298.15 with units K'''
Pr = 1
'''The standard state pressure of 1 bar'''
E_PrTr = 78.47
'''Epsilon_{P_rT_r} is a unitless constant with value of 78.47'''
bigR = 1.9858775
'''The gas constant with value 1.9858775 cal mol^-1 k^-1'''
Psy = 2600
'''The value of this constant is 2600 bar'''
Theta = 228
'''The value of this temperature is 228 Kelvin'''
Upsilon = -5.79865E-05
'''The value of this constant is -5.79865E-05 K^-1'''
Conversion = 41.8393
'''A conversion factor present in DEW publications'''


# ## Importing the Aqueous Species Table from the Sheet

# In[4]:


[nameLst, symbolDict, delGf, delHf, entropy, volume, specHeat, a1x10, 
 a2x10_2, a3, a4x10_4, c1, c2x10_4, omegax10_5, Z, comments] = json.load(open(aqPath))


# #### Code for adding additional aqueous species

# In[5]:


# nameLst.append('ALANINE,AQ')
# symbolDict['ALANINE,AQ'] = 'C3H7NO2'
# delGf['ALANINE,AQ'] = -88810  
# delHf['ALANINE,AQ'] = -132500
# entropy['ALANINE,AQ'] = 38.83
# volume['ALANINE,AQ'] = 60.4
# specHeat['ALANINE,AQ'] = 33.6
# a1x10['ALANINE,AQ'] =  14.9
# a2x10_2['ALANINE,AQ'] = 1.74
# a3['ALANINE,AQ'] = 7.16
#a4x10_4['ALANINE,AQ'] = -3.69
# c1['ALANINE,AQ'] = 49.5
# c2x10_4['ALANINE,AQ'] = -7
# omegax10_5['ALANINE,AQ'] = 0.18
# Z['ALANINE,AQ'] = 0


# In[6]:


# d =[nameLst, symbolDict, delGf, delHf, entropy, volume, specHeat, a1x10, 
#  a2x10_2, a3, a4x10_4, c1, c2x10_4, omegax10_5, Z, comments]
# json.dump(d, open("aqueousLst.txt",'w'))


# ## Importing the Gas Table from the Sheet

#  

# In[7]:


[GasLst,GasSymb,GasDelGf,GasDelHf,GasEntropy,GasCp,GasA,GasBx103,GasCx10_5, GasT] = json.load(open(gasPath))


# # An Object Class that Can Calculate and Return Parameters for Different Options of the Deep Earth Water Model

# In[8]:
def search(string):
    for item in nameLst:
        if str.lower(string) in str.lower(item):
            print(item)
    for item in GasLst:
        if str.lower(string) in str.lower(item):
            print(item)
    for key in mineralDictionary:
        if str.lower(string) in str.lower(key):
            print(key)


class DEW(object):
    def __init__(self):
       
        # User Option Parameters
        self.ptInput = 'Psat'
        '''The temperature and pressure input, options are Regular, Psat, or custom. Default is regular'''
        
        self.RhoOfWater = 'Z&D 2005'
        '''The density of water equation input, can be Zheng and Duan 2005, Zheng and Duan 2009, or custom. Default is Z&D 2005'''
        
        self.forceCustom = False
        '''The option to force custom Rho for P< 1 kb. Default is False'''
        
        self.dielectricEq = 'Supcrt'
        '''The dielectric equation input. The default is Sverjensky.'''
        
        self.ForceSupcrt = True
        '''The option to force supcrt for P < 5 kb. Default is set to true'''
        self.WaterFreeEq = 'D&H 1978'
        '''The option for the Water free energy equation. Options are D&H 1978, integral, and custom
        Default is Delaney and Hegelson 1978.'''
        self.DisplayVolOpt = True
        '''The option to display volume, default set to true'''
        self.PsatDisplayVol = True
        '''The option to display volume under Psat conditions. Default is set to true.'''
        self.DisplayVol = True
        '''Another display volume option. Default to true.'''
        self.equation = 1
        '''A variable that stores the number of the density of water equation. Needs to be renamed'''
        self.diaEq = 1
        '''A variable that stores the number of dielectric constant equation.'''
        self.psat = True
        '''A variable that stores the Psat option defined by input'''
        self.waterDensity = 1
        '''A variable that stores the number of the density of water equation.'''

        
        # Input Arrays
        self.aqueousInputs = []
        '''The array of aqueous inputs and multipliers defined by a user'''
        self.mineralInputs = []
        '''The array of mineral inputs and multipliers defined by a user'''
        self.gasInputs = []
        '''The array of gas inputs and multipliers defined by a user'''
        self.waterInputs = []
        '''An array that defines if water is used in the input and hOw mUcH wAtEr?'''
        
        # Input Matrices
        self.inGasMat = []
        '''A matrix that stores in gasseous inputs with their properties from the dicitonary inputs'''
        self.inAqMat = []
        '''A matrix that stores in aqueous inputs with their properties from the dicitonary inputs'''
        
        # Output Arrays
        self.aqueousOutputs = []
        '''The array of aqueous outputs and multipliers defined by a user'''
        self.mineralOutputs = []
        '''The array of mineral outputs and multipliers defined by a user'''
        self.gasOutputs = []
        '''The array of gas outputs and multipliers defined by a user'''
        self.waterOutputs = []
        '''An array that defines if water is used in the outputand hOw mUcH wAtEr?'''
        
        # Output Matrices
        self.outGasMat = []
        '''A matrix that stores in gasseous outputs with their properties from the dicitonary inputs'''
        self.outAqMat = []
        '''A matrix that stores in aqueous outputs with their properties from the dicitonary inputs'''
        
        # Arrays used for Calculations
        self.tempUsed = []
        '''An array set by the set_TPRho method that contains all the temperatures used for calculation in celsius'''
        self.pressureUsed = []
        '''An array set by the set_TPRho method that contains all the pressures used for calculation'''
        self.tKelvin = []
        '''An array set by the set_TPRho method that contains all the temperatures used for calculation in Kelvin'''
        self.RhoWatArr = []
        '''An array set by the set_TPRho method that contains calculated water densities at the temperatures and pressures used
        '''
        self.DiaArr = []
        '''An array set by the set_TPRho method that contains calculated dielectric constants at temp/pressure used'''
        self.QArr = []
        '''An array set by the set_TPRho method that contains calculated Q constants at temp/pressure used'''
        self.GibbsH2O = []
        '''A collection of the gibbs of water values.'''
        
        # Collections of Custom Values
        self.dielectricCollection = []
        '''If custom values are used for the dielectric constant this will store them to be queried by the custom function'''
        self.gibbsCollection = []
        '''If custom values are used for the gibbs of water this will store them to be queried by the custom function'''
        self.densityCollection = []
        '''If custom values are used for the density of water this will store them to be queried by the custom function'''
        
        # Calculated Matrices
        self.gasInpGibbs = []
        '''Used for debugging, stores the free energy changes of gases'''
        self.aqInpGibbs = []
        '''Used for debugging, stores the free energy changes of aqueous inputs'''
        self.gasInpV = []
        '''Used for debugging, stores the volume changes of gasseous inputs'''
        self.aqInpV = []
        '''Used for debugging, stores the volume changes of aqueous inputs'''
        self.gasOutGibbs = []
        '''Used for debugging, stores the free energy changes of gasseous inputs'''
        self.aqOutGibbs = []
        '''Used for debugging, stores the free energy changes of aqueous outputs'''
        self.gasOutV = []
        '''Used for debugging, stores the volume changes of gasseous outputs'''
        self.aqOutV = []
        '''Used for debugging, stores the volume changes of aqueous outputs'''
        
        #Mineral Matrices
        self.mineralInpGibbs = []
        '''Used for debugging, stores the free energy changes of mineral inputs'''
        self.mineralOutGibbs = []
        '''Used for debugging, stores the free energy changes of mineral outputs'''
        self.mineralInpV = [] 
        '''Used for debugging, stores the volume changes of mineral inputs'''
        self.mineralOutV = []
        '''Used for debugging, stores the volume changes of mineral outputs'''
        
        #Water
        self.waterInpGibbs = []
        '''Used for debugging, stores the free energy changes of water outputs'''
        self.waterInpV = []
        '''Used for debugging, stores the volume changes of water inputs'''
        self.waterOutGibbs = []
        '''Used for debugging, stores the free energy changes of water outputs'''
        self.waterOutV = []
        '''Used for debugging, stores the volume changes of water outputs'''
        
        # Finals Arrays
        self.gibbsLst = []
        '''A storage variable that lists the gibbs free energy changes. Not sure if necessary'''
        self.logK = []
        '''Stores the list of all logK values with temperatures and pressures'''
        self.vLst = []
        '''A storage variable that lists all the volume changes. Not sure if necessary '''
        self.delG = []
        '''Stores the list of all delG values with temperatures and pressures'''
        self.delV = []
        '''Stores the list of all delV values with temperatures and pressures'''
        
        
        # Variables to Help with Plotting
        self.pressRed = []
        '''Reduced pressure list with no repeats'''
        self.tempRed = []
        '''Reduced temperature list with no repeats'''

        self.pLogK = []
        '''LogK split into arrays with respect to the number of isobars'''
        self.pDelG = []
        '''DelG split into arrays with respect to the number of isobars'''
        self.pDelV = []
        '''DelV split into arrays with respect to the number of isobars'''
        
        self.tLogK = []
        '''LogK split into arrays with respect to the number of isotherms'''
        self.tDelG = []
        '''DelG split into arrays with respect to the number of isotherms'''
        self.tDelV = []
        '''DelV split into arrays with respect to the number of isotherms'''
        
        # Variables to run SUPCRTBL
        self.proc = None
        '''Needed to run supcrt'''
        self.pout = None
        '''Needed to run supcrt'''
        self.pin = None
        '''Needed to run supcrt'''
        self.supcrtFile = None
        '''Stores the most recently run SUPCRT file, or none if none have been run'''
        self.supcrtOut = None
        '''Stores the output from calculate_supcrt'''
        

    def clear(self):
        '''Clears variables'''
        self.__init__()
        return
        
    def set_inputs(self):
        '''Call this to set the input Arrays. This is not dependent on anything else being called first.'''
        # A list of integers
        intLst = ['1','2','3','4', '5', '6','7', '8', '9', '10', '11']
        
        # Mineral Loop
        mineralCount = 0
        aqCount = 0
        gasCount = 0
        self.mineralInputs = []
        self.aqueousInputs = []
        self.gasInputs = []
        self.waterInputs = []
        
        while mineralCount < 15:
            mineralCount += 1
            validBool = False
            while not validBool:
                inp = input('Input Mineral Species')
                # can insert mineral validation here if possible
                if inp in mineralDictionary:
                    validBool = True
                elif inp == "":
                    validBool = True
                else:
                    print('Your Species is not in the list, please check your spelling')
                    continue
                validBool2 = False
                
                while not validBool2:
                    inp2 = input('Input Mineral Species Multiplier')
                    if inp2 in intLst:
                        validBool2 = True
                    elif inp == "":
                        validBool2 = True
                    else:
                        print('Your multiplier is invalid, please check to make sure this is an integer')
            if inp == "":
                break
            self.mineralInputs.append([inp, inp2])
            
            
        while aqCount <15:
            aqCount += 1
            validBool = False
            while not validBool:
                inp = input('Input Aqueous Species') 
                if inp in nameLst:
                    validBool = True
                elif inp == "":
                    validBool = True
                else:
                    print('Your Species is not in the list, please check your spelling')
                    continue
                validBool2 = False
                if validBool:
                    while not validBool2:
                        inp2 = input('Input Aqueous Species Multiplier')
                        if inp2 in intLst:
                            validBool2 = True
                        elif inp == "":
                            validBool2 = True
                        else:
                            print('Your multiplier is invalid, please check to make sure this is an integer')
            if inp == "":
                break
            self.aqueousInputs.append([inp, inp2])
            
            
        while gasCount < 15:
            gasCount += 1
            validBool = False
            while not validBool:
                inp = input('Input Gas Species') 
                if inp in GasLst:
                    validBool = True
                elif inp == "":
                    validBool = True
                else:
                    print('Your Species is not in the list, please check your spelling')
                    continue
                if validBool:
                    validBool2 = False
                    while not validBool2:
                        inp2 = input('Input Gas Species Multiplier')
                        if inp2 in intLst:
                            validBool2 = True
                        elif inp == "":
                            validBool2 = True
                        else:
                            print('Your multiplier is invalid, please check to make sure this is an integer')
            if inp == "":
                break
            self.gasInputs.append([inp, inp2])
            
            
            
            # Water
        validBool3 = False
        while not validBool3:
            inpWater = input('Would you like to use water? (yes/no)')
            if inpWater in ['yes', 'no']:
                validBool3 = True
            else:
                print('Please answer yes or no')
                continue
            if inpWater == 'yes':
                validBool3 = False
                while not validBool3:
                    m3 = input('Enter enter water Multiplier')
                    if m3 in intLst:
                        validBool3 = True
                    else:
                        print('Please enter a valid integer multiplier ')
            else: 
                m3 = 0
            self.waterInputs.append([inpWater, m3])
        return
    
    def set_outputs(self):
        '''Call this to set the output Arrays. This is not dependent on anything else being called first.'''
        # A list of integers
        intLst = ['1','2','3','4', '5', '6','7', '8', '9', '10', '11']
        
        # Mineral Loop
        mineralCount = 0
        aqCount = 0
        gasCount = 0
        self.mineralOutputs = []
        self.aqueousOutputs = []
        self.gasOutputs = []
        self.waterOutputs = []


        while mineralCount < 15:
            mineralCount += 1
            validBool = False
            while not validBool:
                inp = input('Output Mineral Species')
                # can insert mineral validation here if possible
    
                validBool = True
        
                validBool2 = False
                while not validBool2:
                    inp2 = input('Output Mineral Species Multiplier')
                    if inp2 in intLst:
                        validBool2 = True
                    elif inp == "":
                        validBool2 = True
                    else:
                        print('Your multiplier is invalid, please check to make sure this is an integer')
            if inp == "":
                break
            self.mineralOutputs.append([inp, inp2])
            
            
        while aqCount <15:
            aqCount += 1
            validBool = False
            while not validBool:
                inp = input('Output Aqueous Species') 
                if inp in nameLst:
                    validBool = True
                elif inp == "":
                    validBool = True
                else:
                    print('Your Species is not in the list, please check your spelling')
                    continue
                validBool2 = False
                if validBool:
                    while not validBool2:
                        inp2 = input('Output Aqueous Species Multiplier')
                        if inp2 in intLst:
                            validBool2 = True
                        elif inp == "":
                            validBool2 = True
                        else:
                            print('Your multiplier is invalid, please check to make sure this is an integer')
            if inp == "":
                break
            self.aqueousOutputs.append([inp, inp2])
            
        while gasCount < 15:
            gasCount += 1
            validBool = False
            while not validBool:
                inp = input('Input Gas Species') 
                if inp in GasLst:
                    validBool = True
                elif inp == "":
                    validBool = True
                else:
                    print('Your Species is not in the list, please check your spelling')
                    continue
                validBool2 = False
                if validBool:
                    while not validBool2:
                        inp2 = input('Input Gas Species Multiplier')
                        if inp2 in intLst:
                            validBool2 = True
                        elif inp == "":
                            validBool2 = True
                        else:
                            print('Your multiplier is invalid, please check to make sure fthis is an integer')
            if inp == "":
                break
            self.gasOutputs.append([inp, inp2])
            
            # Water
        validBool3 = False
        while not validBool3:
            outWater = input('Would you like to use water in the output? (yes/no)')
            if outWater in ['yes', 'no']:
                validBool3 = True
            else:
                print('Please answer yes or no')
            if outWater == 'yes':
                validBool3 = False
                while not validBool3:
                    m3 = input('Enter enter water Multiplier')
                    if m3 in intLst:
                        validBool3 = True
                    else:
                        print('Please enter a valid integer multiplier ')
            else: 
                m3 = 0
            self.waterOutputs.append([outWater, m3])
        return
        
    
    def set_preferences(self):
        '''A function that prompts for user inputs. This is not dependent on anything else being called first. Defaults
        are set to be identical to the example calculation on the Deep Earth Water Model Excel Sheet.'''
        validBool = False
        while not validBool:  
            ptInp = input('Which P-T input would you like to use? "Custom", "Regular", or "Psat"')
            if ptInp in ['Custom', 'Regular', 'Psat']:
                validBool = True
                self.ptInput = ptInp
            else:
                print('Please enter one of the provided options')
       
        validBool = False
        while not validBool:
            RhoOfwater = input('Which density of water would you like to use? "Z&D 2005", "Z&D 2009", or "Custom"')
            if RhoOfwater in ['Z&D 2005', 'Z&D 2009', 'Custom']:
                validBool = True
                self.RhoOfWater = RhoOfwater
            else:
                print('Please enter one of the provided options')
        
        validBool = False
        while not validBool:
            force = input('Force Custom? (yes/no)')
            if force == 'yes':
                validBool = True
            elif force == 'no':
                validBool = True
                self.forceCustom = False
            else:
                print('Please enter one of the provided options')
            
        validBool = False
        while not validBool:
            dia = input('Dielectric Constant Equation Option: "Supcrt", "Franck", "Fernandez", "Sverjensky", or "Custom"')
            if dia in ['Supcrt', 'Franck', 'Fernandez', 'Sverjensky','Custom']:
                validBool = True
                self.dielectricEq = dia
            else:
                print('Please enter one of the provided options')
        
        validBool = False
        while not validBool:
            forceS = input('Force Supcrt? (yes/no)')
            if forceS == 'yes':
                validBool = True
            elif forceS == 'no':
                validBool = True
                self.ForceSupcrt = False
            else:
                print('Please enter one of the provided options')
        
        validBool = False
        while not validBool:
            freeE = input('Water Free Energy Equation Option: "D&H 1978", "Integral", "Custom"')
            if freeE in ['D&H 1978', 'Integral', 'Custom']:
                validBool = True
                self.WaterFreeEq = freeE

        validBool = False
        while not validBool:
            dispO = input('Display Volume Option? (yes/no)')
            if dispO == 'yes':
                validBool = True
            elif dispO == 'no':
                validBool = True
                self.DisplayVolOpt = False
            else:
                print('Please enter one of the provided options')
                 
        validBool = False            
        while not validBool:
            PsatdispO = input('Psat Display Volume Option? (yes/no)')
            if PsatdispO == 'yes':
                validBool = True
            elif PsatdispO == 'no':
                validBool = True
                self.PsatDisplayVol = False
            else:
                print('Please enter one of the provided options')
        
        validBool = False
        while not validBool:
            dispV = input('Display Volume? (yes/no)')
            if dispV == 'yes':
                validBool = True
            elif dispV == 'no':
                validBool = True
                self.DisplayVol = False
            else:
                print('Please enter one of the provided options')
        if self.WaterFreeEq == "Custom" or self.dielectricEq == "Custom" or self.RhoOfWater == "Custom":
            self.dielectricCollection, self.densityCollection, self.gibbsCollection = self.import_custom_sheets()
        return
    
    
    
    
    def import_custom_sheets(self):
        '''A helper function to import custom data from the Deep Earth Water Model.
        This only currently works for an unmodified Deep Earth Water Model Sheet format (6_23_20). 
        This is not dependent on anything else being called first.'''
        
        diaL = pd.read_csv(diePath, header = None)
        dia = diaL.to_numpy()
        dia = dia[4:, 1:]
        diaTrim = dia[1:, 1:]
        diaCollection = []
        for row in range(len(diaTrim)):
            for pressure in range(len(diaTrim[0])):
                # in form pressure, temperature, value
                diaCollection.append([dia[0][pressure + 1], dia[row + 1][0], diaTrim[row][pressure]])

        watDen = pd.read_csv(denPath, header = None)
        w = watDen.to_numpy()
        w = w[4:, 1:]
        wTrim = w[1:,1:]
        watDenCollection = []
        for row in range(len(wTrim)):
            for pressure in range(len(wTrim[0])):
                # in form pressure, temperature, value
                watDenCollection.append([w[0][pressure + 1], w[row + 1][0], wTrim[row][pressure]])

        gibbsOfWater = pd.read_csv(gPath, header = None)
        gibbs = gibbsOfWater.to_numpy()
        gibbs = gibbs[3:,:]
        gibbsTrim = gibbs[1:, 1:]
        gibbsCollection = []
        for row in range(len(gibbsTrim)):
            for pressure in range(len(gibbsTrim[0])):
                # in form pressure, temperature, value
                gibbsCollection.append([gibbs[0][pressure + 1], gibbs[row + 1][0], gibbsTrim[row][pressure]])
        return diaCollection, watDenCollection, gibbsCollection
    
    
    
    def set_TPRho(self):
        '''Sets arrays of temperature, pressure, water density, and Q to be used in the model based on user input. 
        Requires that the input and output arrays have been set up otherwise it will return a divide by 0 error in the 
        calculations.'''
        pressArr = []
        tempArr = []
        self.RhoWatArr = []
        self.DiaArr = []
        self.QArr =[]
        self.gibbsLst = []
        self.logK = []
        self.vLst = []
        self.delG = []
        self.delV = []

        
        if self.ptInput == "Custom":
            ptSheet = pd.read_csv(inpPath,encoding= 'unicode_escape', header = None)
            ptFinder = ptSheet.to_numpy()
            tempArr = [float(i[1]) for i in ptFinder[4:]]
            pressArr = [float(i[0]) for i in ptFinder[4:]]

        elif self.ptInput == "Regular":
            validBool = False
            while not validBool:
                try:
                    templow = int(input('Input the minimum temperature'))
                    temphigh = int(input('Input the maximum temperature'))
                    tempstep = int(input('Input the temperature step'))
                    pmin = float(input('Input the minimum pressure (Kb)'))
                    pmax = float(input('Input the maximum pressure (Kb)'))
                    pstep = float(input('Input the pressure step (Kb)'))
                    validBool = True
                except ValueError:
                    print('You have entered a non-integer value, please start again')
            tempArr = np.arange(start= templow, stop = temphigh + .00001, step = tempstep)
            parrHelp = np.arange(start= pmin, stop = pmax + .00001, step = pstep)
            for i in range(len(parrHelp)):
                pressArr.append([parrHelp[i]]* len(tempArr))
            pressArr = np.multiply(pressArr, 1000)
            tempArr = [tempArr] * len(parrHelp)
            
        elif self.ptInput == "Psat":
            validBool = False
            while not validBool:
                try:
                    templow = int(input('Input the minimum temperature'))
                    temphigh = int(input('Input the mamximum temperature'))
                    tempstep = int(input('Input the temperature step'))
                    validBool = True
                except ValueError:
                    print('You have entered a non-integer value, please start again')
                    
            tempArr = np.arange(start= templow, stop = temphigh + 1, step = tempstep)
            for i in range(len(tempArr)):
                
                if tempArr[i] < 100:
                    pressArr.append(1)
                else:
                    pressArr.append(2.1650906415E-11*np.double(tempArr[i])**5 + 0.0008467019353*np.double(tempArr[i])**2 - 0.17973651666*tempArr[i] + 10.7768850763807)
                
        else:
            # If I've done the checking correctly above it should never reach this
            raise ValueError("You have not set your options yet, please set them before continuing")
        self.tempUsed = np.ndarray.flatten(np.asarray(tempArr))
        self.pressureUsed = np.ndarray.flatten(np.asarray(pressArr))
        self.tKelvin = np.add(self.tempUsed, 273.15)
        
        # code to set options in a way the equations can understand
        if self.ptInput == "Psat":
            self.psat = True
        else:
            self.psat = False
            
        if self.RhoOfWater =='Z&D 2005':
            self.equation = 1
        elif self.RhoOfWater == 'Z&D 2009':
            self.equation = 2
        else:
            self.equation = 3
            
        if self.dielectricEq == "Supcrt":
            self.diaEq = 1
        elif self.dielectricEq == "Franck":
            self.diaEq = 2
        elif self.dielectricEq == "Fernandez":
            self.diaEq = 3
        elif self.dielectricEq == "Sverjensky":
            self.diaEq = 4
        else:
            self.diaEq = 5
        
        # write code to take in custom Rho, G, and Water Values here
        self.densityCollection = np.asarray(self.densityCollection).astype(float)
        self.dielectricCollection = np.asarray(self.dielectricCollection).astype(float)
        self.gibbsCollection = np.asarray(self.gibbsCollection).astype(float)
        
        # Sets the water density array
        for i in range(len(self.pressureUsed)):        
            # For the custom array
            if self.RhoOfWater =="Custom" or (self.forceCustom == True and self.pressureUsed[i] < 1000):
                idx = np.intersect1d(np.where(np.asarray(self.densityCollection) == self.pressureUsed[i]/1000), np.where(np.asarray(self.densityCollection) == self.tempUsed[i]))[0]
                if not np.isnan(self.densityCollection[idx][2]):
                    self.RhoWatArr.append(self.densityCollection[idx][2])
                else:
                    self.RhoWatArr.append(0)
            else:
                self.RhoWatArr.append(DEWEquations.DEWEquations.calculateDensity(self.pressureUsed[i], self.tempUsed[i], self.equation, 0.01, self.psat))
               
        # Sets the dielectric constant array
        for i in range(len(self.pressureUsed)):
            
            # for the custom array
            if self.dielectricEq == "Custom":
                idx = np.intersect1d(np.where(np.asarray(self.dielectricCollection) == self.pressureUsed[i]/1000), np.where(np.asarray(self.dielectricCollection) == self.tempUsed[i]))[0]
                if not np.isnan(self.dielectricCollection[idx][2]):
                    self.DiaArr.append(self.dielectricCollection[idx][2])
                else:
                    self.DiaArr.append(0)
            else:
                if self.ForceSupcrt == True and self.pressureUsed[i] < 5000 and self.psat == False:
                    self.DiaArr.append(DEWEquations.DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], 1, self.psat))
                else:
                    self.DiaArr.append(DEWEquations.DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], self.diaEq, self.psat))
        
        
        ### The function works up until this point, I haven't debugged further yet (6_29_20) ###
        
        # Sets up the Q array
        for i in range(len(self.pressureUsed)):
            if self.DisplayVol == True:
                try:
                    # Has issues with some Q, not sure if problematic
                    self.QArr.append(float(DEWEquations.DEWEquations.calculateQ(self.pressureUsed[i], self.tempUsed[i], self.RhoWatArr[i], self.equation, self.diaEq, self.psat))*np.double(10)**6)
                except:
                    self.QArr.append(0)
            else:
                self.QArr.append(0)
                
        # Sets up custom Gibbs of Water Array:
        if self.WaterFreeEq == "Custom":
            for i in range(len(self.pressureUsed)):
                idx = np.intersect1d(np.where(np.asarray(self.gibbsCollection) == self.pressureUsed[i]/1000), np.where(np.asarray(self.gibbsCollection) == self.tempUsed[i]))[0]
                if not np.isnan(self.gibbsCollection[idx][2]):
                    self.GibbsH2O.append(self.gibbsCollection[idx][2])
                else:
                    self.GibbsH2O.append(0)
        return


    def calculate_matrices(self):
        '''A helper function to aggregate the values to the input and output matrices. 
        It requires both the input and output arrays to be set up to function. It is called within "calculate"'''
        
        self.inAqMat = []
        self.inGasMat = []
        self.outAqMat = []
        self.outGasMat = []
        for i in self.aqueousInputs:
            self.inAqMat.append([i[0],symbolDict[i[0]], delGf[i[0]], delHf[i[0]], entropy[i[0]],volume[i[0]],specHeat[i[0]],
                            a1x10[i[0]], a2x10_2[i[0]], a3[i[0]],a4x10_4[i[0]],c1[i[0]],c2x10_4[i[0]],omegax10_5[i[0]],Z[i[0]], i[1]])
            
        for i in self.gasInputs:
            self.inGasMat.append([i[0],GasSymb[i[0]],GasDelGf[i[0]],GasDelHf[i[0]],GasEntropy[i[0]],GasCp[i[0]], GasA[i[0]],
                             GasBx103[i[0]],GasCx10_5[i[0]],GasT[i[0]], i[1]])
            
        for i in self.aqueousOutputs:
            self.outAqMat.append([i[0],symbolDict[i[0]], delGf[i[0]], delHf[i[0]], entropy[i[0]],volume[i[0]],specHeat[i[0]],
                            a1x10[i[0]], a2x10_2[i[0]], a3[i[0]],a4x10_4[i[0]],c1[i[0]],c2x10_4[i[0]],omegax10_5[i[0]],Z[i[0]], i[1]])

            
        for i in self.gasOutputs:
            self.outGasMat.append([i[0],GasSymb[i[0]],GasDelGf[i[0]],GasDelHf[i[0]],GasEntropy[i[0]],GasCp[i[0]], GasA[i[0]],
                             GasBx103[i[0]],GasCx10_5[i[0]],GasT[i[0]],i[1]])
        return 
    
    def calculate_gas(self):
        '''A helper function to calculate the gasseous columns and output them as a matrix. Specifically returns the arrays 
        gasInGibbs, gasOutGibbs, gasInV, gasOuV. Needs self.tempUsed and self.tKelvin to be set, as well as the input gas matrix.
        It is called within the calculate function.'''
        gasInGibbs = []
        gasOuGibbs = []
        gasInV = []
        gasOuV = []
        for gas in self.inGasMat:
            storelst = []
            storelst2 =[]
            storelst.append(gas[0])
            storelst.append(gas[10])
            storelst2.append(gas[0])
            storelst2.append(gas[10])
            
            for i in range(len(self.tempUsed)):
                if self.DisplayVol == False or self.tempUsed[i] == 0:
                    storelst2.append(0)
                else:
                    storelst2.append(24.465)
                    
            for i in range(len(self.tKelvin)):
                storelst.append(gas[2] - gas[4]*(self.tKelvin[i]-T_r) +                                 gas[6]*(self.tKelvin[i]-T_r - self.tKelvin[i]*np.log(self.tKelvin[i]/T_r)) +                                 gas[7]*(0.001)/2*(2*self.tKelvin[i]*T_r -np.double(self.tKelvin[i])**2 - np.double(T_r) **2) +                                 gas[8]*100000*(np.double(self.tKelvin[i])**2 + np.double(T_r)**2 -2*self.tKelvin[i]*T_r)/(2*self.tKelvin[i]*np.double(T_r)**2))
            gasInGibbs.append(storelst)
            gasInV.append(storelst2)
            
        for gas in self.outGasMat:
            storelst = []
            storelst2 = []
            
            storelst.append(gas[0])
            storelst.append(gas[10])
            storelst2.append(gas[0])
            storelst2.append(gas[10])
            
            for i in range(len(self.tempUsed)):
                if self.DisplayVol == False or self.tempUsed[i] == 0:
                    storelst2.append(0)
                else:
                    storelst2.append(24.465)
                    
            for i in range(len(self.tKelvin)):
                storelst.append(gas[2] - gas[4]*(self.tKelvin[i]-T_r) +                                 gas[6]*(self.tKelvin[i]-T_r - self.tKelvin[i]*np.log(self.tKelvin[i]/T_r)) +                                 gas[7]*(0.001)/2*(2*self.tKelvin[i]*T_r -np.double(self.tKelvin[i])**2 - np.double(T_r) **2) +                                 gas[8]*100000*(np.double(self.tKelvin[i])**2 + np.double(T_r)**2 -2*self.tKelvin[i]*T_r)/(2*self.tKelvin[i]*np.double(T_r)**2))
            gasOuGibbs.append(storelst)
            gasOuV.append(storelst2)
        if len(gasInGibbs) == 0:
            gasInGibbs = [np.zeros(len(self.tKelvin) + 2)]
        if len(gasOuGibbs) == 0:
            gasOuGibbs = [np.zeros(len(self.tKelvin) + 2)]
        if len(gasInV) == 0:
            gasInV = [np.zeros(len(self.tKelvin) + 2)]
        if len(gasOuV) == 0:
            gasOuV = [np.zeros(len(self.tKelvin) + 2)]
        return gasInGibbs, gasOuGibbs, gasInV, gasOuV
    

    
    def calculate_H2O(self):
        '''This function requires input and output matrices to be set. This is called within the calculate function.'''
        waterMatInGibbs = []
        waterMatOutGibbs = []
        waterMatInV = []
        waterMatOutV = []
        if self.WaterFreeEq == 'D&H 1978':
            self.waterDensity = 1
        elif self.WaterFreeEq == 'Integral':
            self.waterDensity = 2
        else:
            self.waterDensity = 3
        
        if self.waterInputs[0][0] == 'yes':
            waterLst = []
            waterLst2 = []
            waterLst.append('H2O')
            waterLst.append(self.waterInputs[0][1])
            waterLst2.append('H2O')
            waterLst2.append(self.waterInputs[0][1])
                                  
            for i in range(len(self.pressureUsed)):
            #for i in range(len(self.pressureUsed)):
                if self.WaterFreeEq == 'Custom':
                    try:
                        if self.GibbsH2O[i] == 0:
                            waterLst.append(0)
                        else:
                            waterLst.append(GibbsH2O[i])
                    except:
                        waterLst.append(GibbsH2O[i])
                else:
                    store = DEWEquations.DEWEquations.calculateGibbsOfWater(self.pressureUsed[i], self.tempUsed[i], self.waterDensity, self.equation, self.psat)
                    waterLst.append(store)
                if self.DisplayVol == True:
                    try:
                        waterLst2.append(18.01528/self.RhoWatArr[i])
                    except:
                        waterLst2.append(0)
                        continue
                else:
                    waterLst2.append(0)
                    
            waterMatInGibbs.append(waterLst)
            waterMatInV.append(waterLst2)
            
        if self.waterOutputs[0][0] =='yes':
            waterLst = []
            waterLst2 = []
            waterLst.append('H2O')
            waterLst.append(self.waterOutputs[0][1])
            waterLst2.append('H2O')
            waterLst2.append(self.waterOutputs[0][1])
            for i in range(len(self.pressureUsed)):
                if self.WaterFreeEq == 'Custom':
                    try:
                        if GibbsH2O[i] == 0:
                            waterLst.append(0)
                        else:
                            waterLst.append(GibbsH2O[i])
                    except:
                        waterLst.append(GibbsH2O[i])
                else:
                    waterLst.append(DEWEquations.DEWEquations.calculateGibbsOfWater(self.pressureUsed[i], self.tempUsed[i], self.waterDensity, self.equation, self.psat))
                if self.DisplayVol == True:
                    try:
                        waterLst2.append(18.01528/self.RhoWatArr[i])
                    except:
                        waterLst2.append(0)
                else:
                    waterLst2.append(0)
                    
            waterMatOutGibbs.append(waterLst)
            waterMatOutV.append(waterLst2)
        if len(waterMatInGibbs) == 0:
            waterMatInGibbs = np.zeros((len(self.tKelvin) + 2))
        if len(waterMatInV) == 0:
            waterMatInV = np.zeros((len(self.tKelvin) + 2))
        if len(waterMatOutGibbs) == 0:
            waterMatOutGibbs = np.zeros((len(self.tKelvin) + 2))
        if len(waterMatOutV) == 0:
            waterMatOutV = np.zeros((len(self.tKelvin) + 2))
            
        return waterMatInGibbs, waterMatInV, waterMatOutGibbs, waterMatOutV
    

    
    def calculate_aq(self):
        '''A helper function to calculate the aqueous columns and output them as a matrix. This is called within calculate.'''
        aqInGibbs = []
        aqOuGibbs = []
        aqInV = []
        aqOuV = []
        for aq in self.inAqMat:
            storelst = []
            storelst2= []
            storelst.append(aq[0])
            storelst.append(aq[15])
            storelst2.append(aq[0])
            storelst2.append(aq[15])
            for i in range(len(self.tKelvin)):
                storelst.append(aq[2] - aq[4] * (self.tKelvin[i] - T_r)
                                - aq[11] * (self.tKelvin[i] * np.log(self.tKelvin[i]/T_r) - self.tKelvin[i] + T_r)
                                - aq[12]*(10**4)*(((1/(self.tKelvin[i]-Theta)) - (1/(T_r-Theta)))*((Theta-self.tKelvin[i])/(Theta))- (self.tKelvin[i]/(Theta*Theta)) * np.log((T_r*(self.tKelvin[i]-Theta))/(self.tKelvin[i]*(T_r-Theta))))
                                + aq[7]*(10**-1)*(self.pressureUsed[i]-Pr)
                                + aq[8]*(10**2)*np.log((Psy+self.pressureUsed[i])/(Psy+Pr))
                                + (1/(self.tKelvin[i]-Theta))*(aq[9]*(self.pressureUsed[i]-Pr)
                                                               + aq[10]*(10**4)*np.log((Psy+self.pressureUsed[i])/(Psy+Pr)))
                                + DEWEquations.DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*((1/self.DiaArr[i])-1)
                                - aq[13]*(10**5)*((1/E_PrTr)-1)
                                + aq[13]*(10**5)*Upsilon*(self.tKelvin[i]-T_r))
                
            for i in range(len(self.pressureUsed)):
                storelst2.append((aq[7]/10 + aq[8]*100/(Psy+self.pressureUsed[i])
                                  + (aq[9] + aq[10]*10000/(Psy+self.pressureUsed[i]))/(self.tKelvin[i]-Theta)
                                  - DEWEquations.DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*(self.QArr[i]*10**-6 )
                                  + (1/self.DiaArr[i] - 1) * DEWEquations.DEWEquations.calculate_domegadP(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14],self.equation,self.psat))*41.84)
                
            aqInGibbs.append(storelst)
            aqInV.append(storelst2)
                                 
        for aq in self.outAqMat:
            storelst = []
            storelst2= []
            storelst.append(aq[0])
            storelst.append(aq[15])
            storelst2.append(aq[0])
            storelst2.append(aq[15])
            for i in range(len(self.tKelvin)):
                storelst.append(aq[2] - aq[4] * (self.tKelvin[i] - T_r)
                                - aq[11] * (self.tKelvin[i] * np.log(self.tKelvin[i]/T_r) - self.tKelvin[i] + T_r)
                                - aq[12]*(10**4)*(((1/(self.tKelvin[i]-Theta)) - (1/(T_r-Theta)))*((Theta-self.tKelvin[i])/(Theta))- (self.tKelvin[i]/(Theta*Theta)) * np.log((T_r*(self.tKelvin[i]-Theta))/(self.tKelvin[i]*(T_r-Theta))))
                                + aq[7]*(10**-1)*(self.pressureUsed[i]-Pr)
                                + aq[8]*(10**2)*np.log((Psy+self.pressureUsed[i])/(Psy+Pr))
                                + (1/(self.tKelvin[i]-Theta))*(aq[9]*(self.pressureUsed[i]-Pr)
                                                               + aq[10]*(10**4)*np.log((Psy+self.pressureUsed[i])/(Psy+Pr)))
                                + DEWEquations.DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*((1/self.DiaArr[i])-1)
                                - aq[13]*(10**5)*((1/E_PrTr)-1)
                                + aq[13]*(10**5)*Upsilon*(self.tKelvin[i]-T_r))
                
            for i in range(len(self.pressureUsed)):
                storelst2.append((aq[7]/10 + aq[8]*100/(Psy+self.pressureUsed[i])
                                  + (aq[9] + aq[10]*10000/(Psy+self.pressureUsed[i]))/(self.tKelvin[i]-Theta)
                                  - DEWEquations.DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*(self.QArr[i]*10**-6 )
                                  + (1/self.DiaArr[i] - 1) * DEWEquations.DEWEquations.calculate_domegadP(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14],self.equation,self.psat))*41.84)
            aqOuGibbs.append(storelst)
            aqOuV.append(storelst2)
        if len(aqInGibbs) == 0:
            aqInGibbs = [np.zeros(len(self.tKelvin) + 2)]
        if len(aqOuGibbs) == 0:
            aqOuGibbs = [np.zeros(len(self.tKelvin) + 2)]
        if len(aqInV) == 0:
            aqInV = [np.zeros(len(self.tKelvin) + 2)]
        if len(aqOuV) == 0:
            aqOuV = [np.zeros(len(self.tKelvin) + 2)]
        return aqInGibbs, aqOuGibbs, aqInV, aqOuV

    
    def calculate(self):
        '''The function called that will update all of the parameters. It has no outputs, but allows certain arrays to be queried.
        '''
        self.calculate_matrices()
        self.waterInpGibbs, self.waterInpV, self.waterOutGibbs, self.waterOutV = self.calculate_H2O()
        self.aqInpGibbs, self.aqOutGibbs, self.aqInpV, self.aqOutV = self.calculate_aq()
        self.gasInpGibbs, self.gasOutGibbs, self.gasInpV, self.gasOutV = self.calculate_gas()
        

        G1 = np.delete(np.asarray(self.waterInpGibbs), [0,1]).astype(np.float) * int(self.waterInputs[0][1])
        V1 = np.delete(np.asarray(self.waterInpV), [0,1]).astype(np.float) * int(self.waterInputs[0][1])
        G4 = np.delete(np.asarray(self.waterOutGibbs), [0,1]).astype(np.float) * int(self.waterOutputs[0][1])
        V4 = np.delete(np.asarray(self.waterOutV), [0,1]).astype(np.float) * int(self.waterOutputs[0][1])
        
        # Gas Loops
        G3, V3 = ([], [])
        for i in range(len(self.gasInpGibbs)):
            G3.append(np.multiply(np.delete(np.asarray(self.gasInpGibbs[i]), [0,1]).astype(np.float), int(self.gasInpGibbs[i][1])))
            V3.append(np.multiply(np.delete(np.asarray(self.gasInpV[i]), [0,1]).astype(np.float), int(self.gasInpV[i][1])))
        G3 = np.sum(G3, axis = 0)
        V3 = np.sum(V3, axis = 0)
        
        G6, V6 = ([], [])
        for i in range(len(self.gasOutGibbs)):
            G6.append(np.multiply(np.delete(np.asarray(self.gasOutGibbs[i]), [0,1]).astype(np.float), int(self.gasOutGibbs[i][1])))
            V6.append(np.multiply(np.delete(np.asarray(self.gasOutV[i]), [0,1]).astype(np.float),  int(self.gasOutV[i][1])))
        G6 = np.sum(G6, axis = 0)
        V6 = np.sum(V6, axis = 0)
        
        # Aqueous Inputs
        G2, V2 = ([], [])
        for i in range(len(self.aqInpGibbs)):
            G2.append(np.multiply(np.delete(np.asarray(self.aqInpGibbs[i]), [0,1]).astype(np.float),  int(self.aqInpGibbs[i][1])))
            V2.append(np.multiply(np.delete(np.asarray(self.aqInpV[i]), [0,1]).astype(np.float),  int(self.aqInpV[i][1])))
        G2 = np.sum(G2, axis = 0)
        V2 = np.sum(V2, axis = 0)    
            
        G5, V5 = ([], [])
        for i in range(len(self.aqOutGibbs)):
            G5.append(np.multiply(np.delete(np.asarray(self.aqOutGibbs[i]), [0,1]).astype(np.float), int(self.aqOutGibbs[i][1])))
            V5.append(np.multiply(np.delete(np.asarray(self.aqOutV[i]), [0,1]).astype(np.float), int(self.aqOutV[i][1])))
        G5 = np.sum(G5, axis = 0)
        V5 = np.sum(V5, axis = 0)

        dG = [np.sum([G4, G5, G6], axis = 0) - np.sum([G1, G2, G3], axis = 0)]
        dV = [np.sum([V4, V5, V6], axis = 0) - np.sum([V1, V2, V3], axis = 0)]
        
        # Adding the mineral contributions if they exist, must be at the same temperatures and pressures 
        if len(self.mineralInputs) > 0:
            for i in range(len(self.mineralInputs)):
                if self.psat == False:
                    myMinPath = mineralDictionary2
                else: 
                    myMinPath = mineralDictionary
                for temp in self.tempUsed:
                    self.mineralInpGibbs.append(np.multiply(myMinPath[self.mineralInputs[i][0]]['delG'][myMinPath[self.mineralInputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
                    self.mineralInpV.append(np.multiply(myMinPath[self.mineralInputs[i][0]]['delV'][myMinPath[self.mineralInputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
            dG = np.sum([dG, np.sum([self.mineralInpGibbs], axis = 0)], axis = 0)
            dV = np.sum([dV, np.sum([self.mineralInpV], axis = 0)], axis = 0)     
            
        if len(self.mineralOutputs) > 0:
            for i in range(len(self.mineralOutputs)):
                for temp in self.tempUsed:
                    self.mineralOutGibbs.append(np.multiply(myMinPath[self.mineralOutputs[i][0]]['delG'][myMinPath[self.mineralOutputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
                    self.mineralOutV.append(np.multiply(myMinPath[self.mineralOutputs[i][0]]['delV'][myMinPath[self.mineralOutputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
            dG = np.sum([dG, -np.sum([self.mineralOutGibbs],axis = 0)], axis = 0)
            dV = np.sum([dV, -np.sum([self.mineralOutV],axis = 0)], axis = 0)  
            
        self.logK = []
        self.delG = []
        self.delV = []
        for i in range(len(dG[0])):
            self.logK.append([-dG[0][i]/(2.302585*self.tKelvin[i]*bigR), self.tempUsed[i], self.pressureUsed[i]])
            self.delG.append([dG[0][i], self.tempUsed[i], self.pressureUsed[i]])
            self.delV.append([dV[0][i], self.tempUsed[i], self.pressureUsed[i]])
            
            
        # Sets plotting arrays for convenient plotting of isotherms/isobars    
        if self.ptInput!= 'Psat' or self.psat == False:
            self.pressRed = list(set(self.pressureUsed))
            self.tempRed = list(set(self.tempUsed))
            self.pressRed.sort()
            self.tempRed.sort()
            
            temppLogK = defaultdict(list)
            temppDelG = defaultdict(list)
            temppDelV = defaultdict(list)
            temptLogK = defaultdict(list)
            temptDelG = defaultdict(list)
            temptDelV = defaultdict(list)

            for logK, temp, pressure in self.logK:
                temppLogK[pressure].append(logK)
                temptLogK[temp].append(logK)

            for delG, temp, pressure in self.delG:
                temppDelG[pressure].append(delG)
                temptDelG[temp].append(delG)

            for delV, temp, pressure in self.delV:
                temppDelV[pressure].append(delV)
                temptDelV[temp].append(delV)

            
            for item in temppDelG:
                self.pDelG.append(temppDelG[item])
            for item in temppDelV:
                self.pDelV.append(temppDelV[item])
            for item in temppLogK:
                self.pLogK.append(temppLogK[item])
                
            for item in temptDelG:
                self.tDelG.append(temptDelG[item])
            for item in temptDelV:
                self.tDelV.append(temptDelV[item])
            for item in temptLogK:
                self.tLogK.append(temptLogK[item])
              
        return

       
       
###############################       
####### Methods to auto #######        
###############################

    def set_tp(self, pt_arr):
        '''Setting the PT values, but automated for the helperfunction. Can also be used to quick set tp with prompted input'''
        pressArr = []
        tempArr = []
        self.RhoWatArr = []
        self.DiaArr = []
        self.QArr =[]
        self.gibbsLst = []
        self.logK = []
        self.vLst = []
        self.delG = []
        self.delV = []
        
        if self.ptInput == "Custom":
            ptSheet = pd.read_excel(DEW_Location, sheet_name = 'Input', header = None)
            ptFinder = ptSheet.to_numpy()
            pressArr = ptFinder[:,79][5:]
            tempArr = ptFinder[:,80][5:]
            storeidx = 0
            storeidxP = 0
            for i in range(len(tempArr)):
                if np.isnan(tempArr[i]) == True:
                    storeidx = int(i)
                    break
            for i in range(len(pressArr)):
                if np.isnan(pressArr[i]) == True:
                    storeidxP = int(i)
                    break

            tempArr = tempArr[:storeidx]
            pressArr = pressArr[:storeidxP]

        elif self.ptInput == "Regular":
            try:
                templow = pt_arr[0][0]
                temphigh =pt_arr[0][1]
                tempstep = pt_arr[0][2]
                pmin = pt_arr[1][0]
                pmax = pt_arr[1][1]
                pstep = pt_arr[1][2]
            except ValueError:
                print('Your PT array is not formatted correctly. Please use the format [[tmin, tmax, tstep][pmin, pmax, pstep]]')
            tempArr = np.arange(start= templow, stop = temphigh + 1, step = tempstep)
            parrHelp = np.arange(start= pmin, stop = pmax + 1, step = pstep)
            for i in range(len(parrHelp)):
                pressArr.append([parrHelp[i]]* len(tempArr))
            pressArr = np.multiply(pressArr, 1000)
            tempArr = [tempArr] * len(parrHelp)
            
        elif self.ptInput == "Psat":
            try:
                templow = pt_arr[0]
                temphigh = pt_arr[1]
                tempstep = pt_arr[2]
                validBool = True
            except ValueError:
                print('Your input is not formatted correctly. Please use the format for psat of [tmin, tmax, tstep]')
                    
            tempArr = np.arange(start= templow, stop = temphigh + 1, step = tempstep)
            for i in range(len(tempArr)):
                
                if tempArr[i] < 100:
                    pressArr.append(1)
                else:
                    pressArr.append(2.1650906415E-11*np.double(tempArr[i])**5 + 0.0008467019353*np.double(tempArr[i])**2 - 0.17973651666*tempArr[i] + 10.7768850763807)
                
        else:
            # If I've done the checking correctly above it should never reach this
            raise ValueError("You have not set your options yet, please set them before continuing")
        self.tempUsed = np.ndarray.flatten(np.asarray(tempArr))
        self.pressureUsed = np.ndarray.flatten(np.asarray(pressArr))
        self.tKelvin = np.add(self.tempUsed, 273.15)
        
        # code to set options in a way the equations can understand
        if self.ptInput == "Psat":
            self.psat = True
        else:
            self.psat = False
            
        if self.RhoOfWater =='Z&D 2005':
            self.equation = 1
        elif self.RhoOfWater == 'Z&D 2009':
            self.equation = 2
        else:
            self.equation = 3
            
        if self.dielectricEq == "Supcrt":
            self.diaEq = 1
        elif self.dielectricEq == "Franck":
            self.diaEq = 2
        elif self.dielectricEq == "Fernandez":
            self.diaEq = 3
        elif self.dielectricEq == "Sverjensky":
            self.diaEq = 4
        else:
            self.diaEq = 5
        
        # write code to take in custom Rho, G, and Water Values here
        self.densityCollection = np.asarray(self.densityCollection).astype(float)
        self.dielectricCollection = np.asarray(self.dielectricCollection).astype(float)
        self.gibbsCollection = np.asarray(self.gibbsCollection).astype(float)
        
        # Sets the water density array
        for i in range(len(self.pressureUsed)):        
            # For the custom array
            if self.RhoOfWater =="Custom" or (self.forceCustom == True and self.pressureUsed[i] < 1000):
                idx = np.intersect1d(np.where(np.asarray(self.densityCollection).astype(float) == self.pressureUsed[i]/1000), np.where(np.asarray(self.densityCollection).astype(float) == self.tempUsed[i]))[0]
                if not np.isnan(self.densityCollection[idx][2]):
                    self.RhoWatArr.append(self.densityCollection[idx][2])
                else:
                    self.RhoWatArr.append(0)
            else:
                self.RhoWatArr.append(DEWEquations.DEWEquations.calculateDensity(self.pressureUsed[i], self.tempUsed[i], self.equation, 0.01, self.psat))
               
        # Sets the dielectric constant array
        for i in range(len(self.pressureUsed)):
            
            # for the custom array
            if self.dielectricEq == "Custom":
                idx = np.intersect1d(np.where(np.asarray(self.dielectricCollection).astype(float) == self.pressureUsed[i]/1000), np.where(np.asarray(self.dielectricCollection).astype(float) == self.tempUsed[i]))[0]
                if not np.isnan(self.dielectricCollection[idx][2]):
                    self.DiaArr.append(self.dielectricCollection[idx][2])
                else:
                    self.DiaArr.append(0)
            else:
                if self.ForceSupcrt == True and self.pressureUsed[i] < 5000 and self.psat == False:
                    self.DiaArr.append(DEWEquations.DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], 1, self.psat))
                else:
                    self.DiaArr.append(DEWEquations.DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], self.diaEq, self.psat))
        
        
        # Sets up the Q array
        for i in range(len(self.pressureUsed)):
            if self.DisplayVol == True:
                try:
                    self.QArr.append(float(DEWEquations.DEWEquations.calculateQ(self.pressureUsed[i], self.tempUsed[i], self.RhoWatArr[i], self.equation, self.diaEq, self.psat))*np.double(10)**6)
                except:
                    self.QArr.append(0)
            else:
                self.QArr.append(0)
                
        # Sets up custom Gibbs of Water Array:
        if self.WaterFreeEq == "Custom":
            for i in range(len(self.pressureUsed)):
                idx = np.intersect1d(np.where(np.asarray(self.gibbsCollection).astype(float) == self.pressureUsed[i]/1000), np.where(np.asarray(self.gibbsCollection).astype(float) == self.tempUsed[i]))[0]
                if not np.isnan(self.gibbsCollection[idx][2]):
                    self.GibbsH2O.append(self.gibbsCollection[idx][2])
                else:
                    self.GibbsH2O.append(0)
        return
    def run(self, pt_arr, min_inp =[], aq_inp = [], g_inp = [], h2o_inp = 0, min_out = [],aq_out =[], g_out = [],h2o_out = 0, 
        ptInp = 'Psat', rhoWat = 'Z&D 2005', forceBool = False, dieEQ = 'Supcrt', forceSC = True, 
        WFEQ ='D&H 1978', dsV = True, pdsV = True, DV = True, EQ = 1, dEQ = 1, pst = True, mWn = 1, makeP = False):

        if h2o_inp > 0:
            self.waterInputs = [['yes',h2o_inp]]
        else:
            self.waterInputs = [['no',0]]
        if h2o_out > 0:
            self.waterOutputs = [['yes',h2o_out]]
        else:
            self.waterOutputs = [['no',0]]

        self.mineralInputs = min_inp
        self.aqueousInputs = aq_inp
        self.gasInputs = g_inp

        self.mineralOutputs = min_out
        self.aqueousOutputs = aq_out
        self.gasOutputs = g_out


        self.ptInput = ptInp
        self.RhoOfWater = rhoWat
        self.forceCustom = forceBool
        self.dielectricEq = dieEQ     
        self.ForceSupcrt = forceSC
        self.WaterFreeEq =  WFEQ
        self.DisplayVolOpt = dsV
        self.PsatDisplayVol = pdsV
        self.DisplayVol = DV
        self.equation = EQ
        self.diaEq = dEQ
        self.psat = pst
        self.waterDensity = mWn

        # to actually run:
        self.set_tp(pt_arr)
        self.calculate()
        if makeP == True:
            self.make_plots()
        return
       
       
 ###### MAKE PLOTS###########
    
    def make_plots(self):
        '''A final function that the user calls to make the plots possible in the Excel spreadsheet. '''
        plt.clf()
        ###### PSAT PLOTS #######
        if self.psat == True or self.ptInput =='Psat':
            plt.figure()
            plt.plot(self.pressureUsed, [i[0] for i in self.logK])
            plt.xlabel('Pressure (bar)')
            plt.ylabel('LogK')
            plt.title('Pressure vs. LogK Psat Curve')
            plt.show()
            
            plt.figure()
            plt.plot(self.pressureUsed, [i[0] for i in self.delG])
            plt.xlabel('Pressure (bar)')
            plt.ylabel('$\Delta$G')
            plt.title('Pressure vs. $\Delta$G Psat Curve')
            plt.show()
            
            plt.figure()
            plt.plot(self.pressureUsed, [i[0] for i in self.delV])
            plt.xlabel('Pressure (bar)')
            plt.ylabel('$\Delta$V')
            plt.title('Pressure vs. $\Delta$V Psat Curve')
            plt.show()
            
            plt.figure()
            plt.plot(self.tempUsed, [i[0] for i in self.logK])
            plt.xlabel('Temperature ($^\circ$ C)')
            plt.ylabel('LogK')
            plt.title('Temperature vs. LogK Psat Curve')
            plt.show()
            
            plt.figure()
            plt.plot(self.tempUsed, [i[0] for i in self.delG])
            plt.xlabel('Temperature ($^\circ$ C)')
            plt.ylabel('$\Delta$G')
            plt.title('Temperature vs. $\Delta$G Psat Curve')
            plt.show()
            
            plt.figure()
            plt.plot(self.tempUsed, [i[0] for i in self.delV])
            plt.xlabel('Temperature ($^\circ$ C)')
            plt.ylabel('$\Delta$V')
            plt.title('Temperature vs. $\Delta$V Psat Curve')
            plt.show()
            
        ####### NON PSAT PLOTS ########    
        else:
            # T Plots
            plt.figure()
            for i in self.pDelG:
                plt.plot(self.tempRed, i)
                plt.legend(self.pressRed,bbox_to_anchor=(1.05, 1), title = 'Pressure (bar)', loc='upper left')
                plt.xlabel('Temperature ($^\circ$C)')
                plt.ylabel('$\Delta$G')
                plt.title('Temperature vs. $\Delta$G')
            plt.show()
                
            plt.figure()
            for i in self.pDelV:
                plt.plot(self.tempRed, i)
                plt.legend(self.pressRed,bbox_to_anchor=(1.05, 1), title = 'Pressure (bar)', loc='upper left')
                plt.xlabel('Temperature ($^\circ$C)')
                plt.ylabel('$\Delta$V')
                plt.title('Temperature vs. $\Delta$V')
            plt.show()
                          
            plt.figure()
            for i in self.pLogK:
                plt.plot(self.tempRed, i)
                plt.legend(self.pressRed,bbox_to_anchor=(1.05, 1), title = 'Pressure (bar)', loc='upper left')
                plt.xlabel('Temperature ($^\circ$C)')
                plt.ylabel('LogK')
                plt.title('Temperature vs. LogK')
            plt.show()
                          
            # P Plots
            plt.figure()              
            for i in self.tDelG:
                plt.plot(self.pressRed, i)
                plt.legend(self.tempRed,bbox_to_anchor=(1.05, 1), title = 'Temperature ($^\circ$C)', loc='upper left')
                plt.xlabel('Pressure (bar)')
                plt.ylabel('$\Delta$G')
                plt.title('Pressure vs. $\Delta$G')
            plt.show()
            
            plt.figure()    
            for i in self.tDelV:
                plt.plot(self.pressRed, i)
                plt.legend(self.tempRed,bbox_to_anchor=(1.05, 1), title = 'Temperature ($^\circ$C)', loc='upper left')
                plt.xlabel('Pressure (bar)')
                plt.ylabel('$\Delta$V')
                plt.title('Pressure vs. $\Delta$V')
            plt.show()
            
            plt.figure()              
            for i in self.tLogK:
                plt.plot(self.pressRed, i)
                plt.legend(self.tempRed,bbox_to_anchor=(1.05, 1), title = 'Temperature ($^\circ$C)', loc='upper left')
                plt.xlabel('Pressure (bar)')
                plt.ylabel('LogK')
                plt.title('Pressure vs. LogK')
            plt.show()
        return
       
       
       
       
       
       
#############################       
######### OTHER #############
#############################       
       
       
    def export_to_csv(self):
        dV = [row[0] for row in self.delV]
        dG = [row[0] for row in self.delG]
        lK = [row[0] for row in self.logK]
        T = [row[1] for row in self.logK]
        P = [row[2] for row in self.logK]
        output_array = np.column_stack([T,P, dV,dG,lK])
        df = pd.DataFrame(output_array)
        df.columns = ['Temperature','Pressure','delV','delG','LogK']
        name = input('Input the name of the CSV file')
        finalName = name + ".csv"
        df.to_csv(finalName, index = False)
        
    def options(self):
        print('Welcome to DEWPython, here are the options you can run:')
        print('1. DEW(): this initializes a Deep Earth Water Model Object')
        print('  -The DEW object requires the set_inputs, set_outputs, set_TPRho, and calculate methods to be run.')
        print('  -You can also utilize the import_custom_sheets method to import custom CSV data')
        print('  -After calculating you can use the make_plots or export_to_csv methods.')
        print('2. run_supcrt: this initializes an inline run of SUPCRTBL')
        print('  -After initializing the SUPCRTBL object, run calculate_supcrt to store the supcrt outputs in arrays')
        print('  -You can also use run_supcrt on a supcrt ouput file that has already been run by adding the optional argument of the file name')
        print('  -After this you can run make_supcrt_plots to plot the supcrt files akin the a DEW file')
   
   
   

   
   
   
   
   
######################################   
####### METHODS FOR SUPCRT ###########
######################################

    def outLoop(self):
        '''A helper function to allow SUPCRTBL to run'''
        running = True
        while(running):
            line = self.pout.readline().decode(sys.stdout.encoding)
            print(line, end='')
            running='\n' in line
        print('Finished')
    
    def run_supcrt(self, version = '96'):
        '''A function that runs the pre-compiled SUPCRTBL found in the file folder'''
        if version != '96':
            sup_path = '/'.join(('resources', 'SUPCRTBL.exe'))
            supPath =  pkg_resources.resource_filename(resource_package, sup_path)
        self.proc = subprocess.Popen(supPath,shell = True, stdout = subprocess.PIPE, stdin = subprocess.PIPE, stderr = subprocess.STDOUT)
        self.pout = self.proc.stdout
        self.pin = self.proc.stdin
        threading.Thread(target=self.outLoop).start()
        while(self.proc.poll() is None):
            var = input('User Input: ')
            if '.txt' in var:
                self.supcrtFile = op.dirname(op.abspath(__file__)) + '\\resources\\' + var
            inp=bytearray(var +'\n', sys.stdin.encoding)
            if(self.proc.poll() is None):
                self.pin.write(inp)
                self.pin.flush()
        return
       
    def calculate_supcrt_special(self, customFile = None):
        '''Calculates the output from either SUPCRTBL/SUPCRT96 at isothermal/isobaric temperatures'''
        returnLst = {}
       
        if customFile != None:
            filename = op.dirname(op.abspath(os.getcwd()))+ '\\' + customFile
        elif len(self.supcrtFile) ==0:
            raise ValueError("You haven't run SUPCRT yet")
        else:
            filename = self.supcrtFile
          
        with open(filename, 'r') as f:
            impor = f.read()
            import_data = impor.replace('\t', ' ')

        split = import_data.split('\n')
        for i in range(len(split)):
            try:
                if 'ISOTHERMS(degC)' in split[i]:
                    finalTemp = " ".join(split[i].split()).split(' ')[5]
                    finalPress = " ".join(split[i+1].split()).split(' ')[6]
                    returnVar = input('Enter reaction title')
            except:
                continue

            if 'STANDARD STATE PROPERTIES OF THE REACTION AT ELEVATED TEMPERATURES AND PRESSURES' in split[i]:
                subLst = []
                temp = []
                pres = []
                DH2 = []
                lgK = []
                dlG = []
                dlH = []
                dlS = []
                dlV = []
                dlCp = []
                subDict = {}
                for item in split[(i+4):]:
                    if len(item) > 0:
                        subLst = (" ".join(item.split())).split(' ')
                        try:
                            float(subLst[0])
                        except:
                            continue
                        try:
                            a = subLst[0]
                            b = subLst[1]
                            c = subLst[2]
                            d = subLst[3]
                            e = subLst[4]
                            f = subLst[5]
                            g = subLst[6]
                            h = subLst[7]
                            i = subLst[8]
                            temp.append(a)
                            pres.append(b)
                            DH2.append(c)
                            lgK.append(d)
                            dlG.append(e)
                            dlH.append(f)
                            dlS.append(g) 
                            dlV.append(h)
                            dlCp.append(i)
                            if float(subLst[0]) == finalTemp and float(subLst[1]) == finalPress:
                                break
                        except:
                            continue

                subDict['Temperature'] = [float(i) for i in temp]
                subDict['Pressure'] = [float(i) for i in pres]
                DH2Lst = []
                lgKLst = []
                dlGLst = []
                dlHLst = []
                dlSLst = []
                dlVLst = []
                dlCpLst = []
                for i in range(len(DH2)):
                    try: 
                        DH2Lst.append(float(DH2[i]))
                    except:
                        DH2Lst.append(0)

                    try:
                        lgKLst.append(float(lgK[i]))
                    except: 
                        lgKLst.append(0)

                    try:
                        dlGLst.append(float(dlG[i]))
                    except:
                        dlGLst.append(0)

                    try:
                        dlHLst.append(float(dlH[i]))
                    except:
                        dlHLst.append(0)

                    try:
                        dlSLst.append(float(dlS[i]))
                    except:
                        dlSLst.append(0)

                    try:
                        dlVLst.append(float(dlV[i]))
                    except:
                        dlVLst.append(0)

                    try:
                        dlCpLst.append(dlCp[i])
                    except:
                        dlCpLst.append[0]

                subDict['DH2O'] = DH2Lst
                subDict['LogK'] = lgKLst
                subDict['delG'] = dlGLst
                subDict['delH'] = dlHLst
                subDict['delS'] = dlSLst
                subDict['delV'] = dlVLst
                subDict['delCp'] = dlCpLst
                returnLst[returnVar] = subDict
            self.supcrtOut = returnLst
    
    def supcrt_inp(self, rxn_lst, reaction_type = 'psat'):
        '''Takes a list of reaction lists (comprised of tuples) and runs supcrt'''
        for reaction in rxn_lst:
            proc = subprocess.Popen('supcrt96.exe',stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)
            pout = proc.stdout
            pin = proc.stdin
            it = 0
            rxnVar = 'realReac.con'

            if reaction_type != 'psat':
                rxnVar = 'Xtend.con'

            title = input('What is the title of your reaction?')
            comm = ['n', 'updateSlop1.dat', '2', rxnVar, '2', '1', title]
            for component in reaction:
                if component[1] not in nameLst:
                    print(str(component[1]) + ' is not in the slop16 database. Please check your spelling and try again. You can use the search function the query the database.')
                else:
                    comm.append(str(component[0]) + ' ' + component[1])

            comm.append('0')
            comm.append('y')
            comm.append('n')
            comm.append(title + '.txt')
            comm.append('1')
            comm.append('1')
            comm.append('empty')

            def outLoop():
                running = True
                while(running):
                    line = pout.readline().decode(sys.stdout.encoding)
                    running='\n' in line

            threading.Thread(target=outLoop).start()
            while(proc.poll() is None): 
                try:
                    inp = comm[it]
                    it += 1
                #     inp = bytearray(input('User Input: ')+'\n',sys.stdin.encoding)
                    if(proc.poll() is None):
                        pin.write(bytearray(inp+'\n',sys.stdin.encoding))
                        pin.flush()
                except:
                    pass
        return
        
    def calculate_supcrt(self, customFile = None):
        '''Calculates an output of thermodynamic properties from a SUPCRTBL output file in the same directory. User must input 
        stopping temperature and pressure to allow the program to calculate properly.
        '''
        returnLst = {}
        max_temp = input('Input the Maximum Temperature')
        max_press = input('Input the Maximum Pressure')
        max_temp = float(max_temp)
        max_press = float(max_press)
        if customFile != None:
            file_Path = op.dirname(op.abspath(os.getcwd()))+ '\\' + customFile
        elif len(self.supcrtFile) ==0:
            raise ValueError("You haven't run SUPCRT yet")
        else:
            filename = self.supcrtFile
            filePath='/'.join(('resources', filename))
            file_Path =  pkg_resources.resource_filename(resource_package, filepath)
        with open(file_Path, 'r') as f:
            impor = f.read()
            import_data = impor.replace('\t', ' ')

        split = import_data.split('\n')
        for i in range(len(split)):
            try:
                if 'REACTION TITLE' in split[i]:
                    returnVar = " ".join(split[i+1].split()).split(' ')[0]
                elif '************************************ REACTION' in split[i]:
                    returnVar = " ".join(split[i+4].split()).split(' ')[0]
            except:
                continue
            if 'STANDARD STATE PROPERTIES OF THE REACTION AT ELEVATED TEMPERATURES AND PRESSURES' in split[i]:
                subLst = []
                temp = []
                pres = []
                DH2 = []
                lgK = []
                dlG = []
                dlH = []
                dlS = []
                dlV = []
                dlCp = []
                subDict = {}
                for item in split[(i+7):]:
                    try:
                        if len(item) > 0:
                            subLst = (" ".join(item.split())).split(' ')
                            temp.append(subLst[0])
                            pres.append(subLst[1])
                            DH2.append(subLst[2])
                            lgK.append(subLst[3])
                            dlG.append(subLst[4])
                            dlH.append(subLst[5])
                            dlS.append(subLst[6]) 
                            dlV.append(subLst[7])
                            dlCp.append(subLst[8])   
                        if float(subLst[0]) == max_temp and float(subLst[1]) == max_press:
                            break
                    except:
                        continue
                subDict['Temperature'] = [float(i) for i in temp]
                subDict['Pressure'] = [float(i) for i in pres]
                subDict['DH2O'] = [float(i) for i in DH2]
                subDict['LogK'] = [float(i) for i in lgK]
                subDict['delG'] = [float(i) for i in dlG]
                subDict['delH'] = [float(i) for i in dlH]
                subDict['delS'] = [float(i) for i in dlS]
                storeLst = []
                for i in dlV:
                    if  i =='*********' or i =='NaN':
                        storeLst.append(0)
                    else:
                        storeLst.append(i)
                subDict['delV'] = storeLst
                subDict['delCp'] = [float(i) for i in dlCp]
                returnLst[returnVar] = subDict
            self.supcrtOut = returnLst
            
    def make_supcrt_plots(self):
        '''Creates plots of LogK and delV for already-calculated SUPCRTBL functions. Produces the same set of plots as the DEW produces'''
        for i in self.supcrtOut:
            plt.figure()
            plt.plot(self.supcrtOut[i]['Temperature'], self.supcrtOut[i]['LogK'])
            plt.title('LogK vs. Temp for ' + i)
            plt.xlabel('Temp, Deg C')
            plt.ylabel('LogK')
            plt.show()
            plt.figure()
            plt.plot(self.supcrtOut[i]['Pressure'], self.supcrtOut[i]['LogK'])
            plt.title('LogK vs. Pressure for ' + i)
            plt.ylabel('LogK')
            plt.xlabel('Pressure (Kb)')
            plt.show()  
            
        for i in self.supcrtOut:
            plt.figure()
            plt.plot(self.supcrtOut[i]['Temperature'], self.supcrtOut[i]['delG'])
            plt.title('delV vs. Temp for ' + i)
            plt.xlabel('Temp, Deg C')
            plt.ylabel('delG')
            plt.show()
            plt.figure()
            plt.plot(self.supcrtOut[i]['Pressure'], self.supcrtOut[i]['delG'])
            plt.title('delV vs. Pressure for ' + i)
            plt.ylabel('delG')
            plt.xlabel('Pressure (Kb)')
            plt.show()
            
        for i in self.supcrtOut:
            plt.figure()
            plt.plot(self.supcrtOut[i]['Temperature'], self.supcrtOut[i]['delV'])
            plt.title('delV vs. Temp for ' + i)
            plt.xlabel('Temp, Deg C')
            plt.ylabel('delV')
            plt.show()
            plt.figure()
            plt.plot(self.supcrtOut[i]['Pressure'], self.supcrtOut[i]['delV'])
            plt.title('delV vs. Temp for ' + i)
            plt.ylabel('delV')
            plt.xlabel('Pressure (Kb)')
            plt.show()


# In[9]:



# In[ ]:




