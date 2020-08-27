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
import os.path as op

# ### Defining a Global Variables (Location and Constants)

# In[3]:
mineralPath = op.dirname(op.abspath(__file__)) + '\\resources\\mineralDictionary.txt'
gasPath = op.dirname(op.abspath(__file__)) + '\\resources\\gasLst.txt'
aqPath = op.dirname(op.abspath(__file__)) + '\\resources\\aqueousLst.txt'
diePath =  op.dirname(op.abspath(__file__)) + '\\resources\\dielectric.csv'
inpPath =  op.dirname(op.abspath(__file__)) + '\\resources\\input.csv'
denPath =  op.dirname(op.abspath(__file__)) + '\\resources\\Wat_den.csv'
gPath =  op.dirname(op.abspath(__file__)) + '\\resources\\water_gibbs.csv'
supPath =  op.dirname(op.abspath(__file__)) + '\\resources\\SUPCRTBL.exe'

global Tr, bigQ, Chi, Pr, E_PrTr, bigR, Psi, Theta, Upsilon, Conversion, mineralDictionary

mineralDictionary = json.load(open(mineralPath))
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
        self.psat = False
        '''A variable that stores the Psat option defined by input'''
        self.myWatNumber = 1
        '''A variable that stores the number of the density of water equation.'''

        
        # Input Arrays
        self.aqueousInputs = []
        '''The array of aqueous inputs and multipliers defined by a user'''
        self.mineralInputs = []
        '''The array of mineral inputs and multipliers defined by a user'''
        self.gasInputs = []
        '''The array of gas inputs and multipliers defined by a user'''
        self.waterInp = []
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
        self.waterOut = []
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
        self.mineralsGInp = []
        '''Used for debugging, stores the free energy changes of mineral inputs'''
        self.mineralsGOutput = []
        '''Used for debugging, stores the free energy changes of mineral outputs'''
        self.mineralsVInp = [] 
        '''Used for debugging, stores the volume changes of mineral inputs'''
        self.mineralsVOutput = []
        '''Used for debugging, stores the volume changes of mineral outputs'''
        
        #Water
        self.InWaterG = []
        '''Used for debugging, stores the free energy changes of water outputs'''
        self.InWaterV = []
        '''Used for debugging, stores the volume changes of water inputs'''
        self.OutWaterG = []
        '''Used for debugging, stores the free energy changes of water outputs'''
        self.OutWaterV = []
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
        
        while mineralCount < 5:
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
            
            
        while aqCount <6:
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
            
            
        while gasCount < 3:
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
        self.inpWater = []
        while not validBool3:
            inpWater = input('Would you like to use water? (yes/no)')
            if inpWater in ['yes', 'no']:
                validBool3 = True
                self.inpWater = inpWater
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
            self.waterInp.append([inpWater, m3])
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
        self.waterOut = []


        while mineralCount < 5:
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
            
            
        while aqCount <6:
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
            
        while gasCount < 3:
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
            self.waterOut.append([outWater, m3])
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
            self.dielectricCollection, self.densityCollection, self.gibbsCollection = import_custom_sheets()
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
        gibbs = gibbs[4:,1:]
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
                    pmin = int(input('Input the minimum pressure'))
                    pmax = int(input('Input the maximum pressure'))
                    pstep = int(input('Input the pressure step'))
                    validBool = True
                except ValueError:
                    print('You have entered a non-integer value, please start again')
            tempArr = np.arange(start= templow, stop = temphigh + 1, step = tempstep)
            parrHelp = np.arange(start= pmin, stop = pmax + 1, step = pstep)
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
        
        # Sets the water density array
        for i in range(len(self.pressureUsed)):        
            # For the custom array
            if self.RhoOfWater =="Custom" or (self.forceCustom == True and self.pressureUsed[i] < 1000):
                idx = np.intersect1d(np.where(np.asarray(self.densityCollection) == pressureUsed[i]/1000), np.where(np.asarray(self.densityCollection) == self.tempUsed[i]))[0]
                if not np.isnan(RhoCollection[idx][2]):
                    self.RhoWatArr.append(self.densityCollection[idx][2])
                else:
                    self.RhoWatArr.append(0)
            else:
                self.RhoWatArr.append(DEWEquations.calculateDensity(self.pressureUsed[i], self.tempUsed[i], self.equation, 0.01, self.psat))
               
        # Sets the dielectric constant array
        for i in range(len(self.pressureUsed)):
            
            # for the custom array
            if self.dielectricEq == "Custom":
                idx = np.intersect1d(np.where(np.asarray(self.dielectricCollection) == pressureUsed[i]/1000), np.where(np.asarray(self.dielectricCollection) == self.tempUsed[i]))[0]
                if not np.isnan(self.dielectricCollection[idx][2]):
                    self.DiaArr.append(self.dielectricCollection[idx][2])
                else:
                    self.DiaArr.append(0)
            else:
                if self.ForceSupcrt == True and self.pressureUsed[i] < 5000 and self.psat == False:
                    self.DiaArr.append(DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], 1, self.psat))
                else:
                    self.DiaArr.append(DEWEquations.calculateEpsilon(self.RhoWatArr[i], self.tempUsed[i], self.diaEq, self.psat))
        
        
        ### The function works up until this point, I haven't debugged further yet (6_29_20) ###
        
        # Sets up the Q array
        for i in range(len(self.pressureUsed)):
            if self.DisplayVol == True:
                try:
                    # Has issues with some Q, not sure if problematic
                    self.QArr.append(float(DEWEquations.calculateQ(self.pressureUsed[i], self.tempUsed[i], self.RhoWatArr[i], self.equation, self.diaEq, self.psat))*np.double(10)**6)
                except:
                    self.QArr.append(0)
            else:
                self.QArr.append(0)
                
        # Sets up custom Gibbs of Water Array:
        if self.WaterFreeEq == "Custom":
            for i in range(len(self.pressureUsed)):
                idx = np.intersect1d(np.where(np.asarray(self.gibbsCollection) == pressureUsed[i]/1000), np.where(np.asarray(self.gibbsCollection) == self.tempUsed[i]))[0]
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
            self.myWatNumber = 1
        elif self.WaterFreeEq == 'Integral of Volume':
            self.myWatNumber = 2
        else:
            self.myWatNumber = 3
        
        if self.waterInp[0][0] == 'yes':
            waterLst = []
            waterLst2 = []
            waterLst.append('H2O')
            waterLst.append(self.waterInp[0][1])
            waterLst2.append('H2O')
            waterLst2.append(self.waterInp[0][1])
                                  
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
                    store = DEWEquations.calculateGibbsOfWater(self.pressureUsed[i], self.tempUsed[i], self.myWatNumber, self.equation, self.psat)
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
            
        if self.waterOut[0][0] =='yes':
            waterLst = []
            waterLst2 = []
            waterLst.append('H2O')
            waterLst.append(self.waterOut[0][1])
            waterLst2.append('H2O')
            waterLst2.append(self.waterOut[0][1])
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
                    waterLst.append(DEWEquations.calculateGibbsOfWater(self.pressureUsed[i], self.tempUsed[i], self.myWatNumber, self.equation, self.psat))
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
                                + DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*((1/self.DiaArr[i])-1)
                                - aq[13]*(10**5)*((1/E_PrTr)-1)
                                + aq[13]*(10**5)*Upsilon*(self.tKelvin[i]-T_r))
                
            for i in range(len(self.pressureUsed)):
                storelst2.append((aq[7]/10 + aq[8]*100/(Psy+self.pressureUsed[i])
                                  + (aq[9] + aq[10]*10000/(Psy+self.pressureUsed[i]))/(self.tKelvin[i]-Theta)
                                  - DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*(self.QArr[i]*10**-6 )
                                  + (1/self.DiaArr[i] - 1) * DEWEquations.calculate_domegadP(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14],self.equation,self.psat))*41.84)
                
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
                                + DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*((1/self.DiaArr[i])-1)
                                - aq[13]*(10**5)*((1/E_PrTr)-1)
                                + aq[13]*(10**5)*Upsilon*(self.tKelvin[i]-T_r))
                
            for i in range(len(self.pressureUsed)):
                storelst2.append((aq[7]/10 + aq[8]*100/(Psy+self.pressureUsed[i])
                                  + (aq[9] + aq[10]*10000/(Psy+self.pressureUsed[i]))/(self.tKelvin[i]-Theta)
                                  - DEWEquations.calculateOmega(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14])*(self.QArr[i]*10**-6 )
                                  + (1/self.DiaArr[i] - 1) * DEWEquations.calculate_domegadP(self.pressureUsed[i],self.tempUsed[i],self.RhoWatArr[i],aq[0],aq[13]*(10**5),aq[14],self.equation,self.psat))*41.84)
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
        self.InWaterG, self.InWaterV, self.OutWaterG, self.OutWaterV = self.calculate_H2O()
        self.aqInpGibbs, self.aqOutGibbs, self.aqInpV, self.aqOutV = self.calculate_aq()
        self.gasInpGibbs, self.gasOutGibbs, self.gasInpV, self.gasOutV = self.calculate_gas()
        

        G1 = np.delete(np.asarray(self.InWaterG), [0,1]).astype(np.float) * int(self.waterInp[0][1])
        V1 = np.delete(np.asarray(self.InWaterV), [0,1]).astype(np.float) * int(self.waterInp[0][1])
        G4 = np.delete(np.asarray(self.OutWaterG), [0,1]).astype(np.float) * int(self.waterOut[0][1])
        V4 = np.delete(np.asarray(self.OutWaterV), [0,1]).astype(np.float) * int(self.waterOut[0][1])
        
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
                for temp in self.tempUsed:
                    self.mineralsGInp.append(np.multiply(mineralDictionary[self.mineralInputs[i][0]]['delG'][mineralDictionary[self.mineralInputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
                    self.mineralsVInp.append(np.multiply(mineralDictionary[self.mineralInputs[i][0]]['delV'][mineralDictionary[self.mineralInputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
            dG = np.sum([dG, np.sum([self.mineralsGInp], axis = 0)], axis = 0)
            dV = np.sum([dV, np.sum([self.mineralsVInp], axis = 0)], axis = 0)     
            
        if len(self.mineralOutputs) > 0:
            for i in range(len(self.mineralOutputs)):
                for temp in self.tempUsed:
                    self.mineralsGOutput.append(np.multiply(mineralDictionary[self.mineralOutputs[i][0]]['delG'][mineralDictionary[self.mineralOutputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
                    self.mineralsVOutput.append(np.multiply(mineralDictionary[self.mineralOutputs[i][0]]['delV'][mineralDictionary[self.mineralOutputs[i][0]]['Temperature'].index(temp)], int(self.mineralInputs[i][1]))) 
            dG = np.sum([dG, -np.sum([self.mineralsGOutput],axis = 0)], axis = 0)
            dV = np.sum([dV, -np.sum([self.mineralsVOutput],axis = 0)], axis = 0)  
            
        self.logK = []
        self.delG = []
        self.delV = []
        for i in range(len(dG[0])):
            self.logK.append([-dG[0][i]/(2.302585*self.tKelvin[i]*bigR), self.tempUsed[i], self.pressureUsed[i]])
            self.delG.append([dG[0][i], self.tempUsed[i], self.pressureUsed[i]])
            self.delV.append([dV[0][i], self.tempUsed[i], self.pressureUsed[i]])
        return
    
    def make_plots(self):
        '''A final function that the user calls to make the plots possible in the DEW Excel spreadsheet. '''
        plt.clf()
        press = list(set(self.pressureUsed))
        temper = list(set(self.tempUsed))
    
        press.sort()
        temper.sort()
        
        pLogK = defaultdict(list)
        pDelG = defaultdict(list)
        pDelV = defaultdict(list)
        tLogK = defaultdict(list)
        tDelG = defaultdict(list)
        tDelV = defaultdict(list)
        
        for logK, temp, pressure in self.logK:
            pLogK[pressure].append(logK)
            tLogK[temp].append(logK)
            
        for delG, temp, pressure in self.delG:
            pDelG[pressure].append(delG)
            tDelG[temp].append(delG)
            
        for delV, temp, pressure in self.delV:
            pDelV[pressure].append(delV)
            tDelV[temp].append(delV)
            
        # Plots for logK
        try:
            pKplot = sorted(pLogK.items()) # sorted by key, return a list of tuples
            x1, y1 = zip(*pKplot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x1, y1)
            if self.psat == False:
                plt.legend(temper, title = "Temperatures (C)")
            plt.xlabel('Pressure (bar)')
            plt.ylabel('LogK')
            plt.title('Pressure vs. LogK')
            plt.show()
        except:
            y1 = list(y1)
            xlst = []
            ylst = []
            for i in range(len(y1)):
                for j in range(len(y1[i])):
                    xlst.append(x1[i])
                    ylst.append(y1[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Pressure (bar)')
            plt.ylabel('LogK')
            plt.title('Pressure vs. LogK Psat Curve')
                    
        plt.figure()
        
        try:
            tKplot = sorted(tLogK.items()) # sorted by key, return a list of tuples
            x2, y2 = zip(*tKplot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x2, y2)
            if self.psat == False:
                plt.legend(press, title = "Pressure (Bar)")
            plt.xlabel('Temperature (C)')
            plt.ylabel('LogK')
            plt.title('Temperature vs. LogK')
            plt.show()
            
        except:
            y2 = list(y2)
            xlst = []
            ylst = []
            for i in range(len(y2)):
                for j in range(len(y2[i])):
                    xlst.append(x2[i])
                    ylst.append(y2[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Temp (C)')
            plt.ylabel('LogK')
            plt.title('Temp vs. LogK Psat Curve')

        plt.figure()
        # Plots for delG
        try:
            pDelGPlot = sorted(pDelG.items()) # sorted by key, return a list of tuples
            x3, y3 = zip(*pDelGPlot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x3, y3)
            if self.psat == False:
                plt.legend(temper, title = "Temperatures (C)")
            plt.xlabel('Pressure (bar)')
            plt.ylabel('Change in Free Energy (DelG)')
            plt.title('Pressure vs. DelG')
            plt.show()
            
        except:
            y3 = list(y3)
            xlst = []
            ylst = []
            for i in range(len(y3)):
                for j in range(len(y3[i])):
                    xlst.append(x3[i])
                    ylst.append(y3[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Pressure (bar)')
            plt.ylabel('DelG')
            plt.title('Pressure vs. DelG Psat Curve')
        
        plt.figure()
        try:
            tDelGPlot = sorted(tDelG.items()) # sorted by key, return a list of tuples
            x4, y4 = zip(*tDelGPlot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x4, y4)
            if self.psat == False:
                plt.legend(press, title = "Pressure (Bar)")
            plt.xlabel('Temperature (C)')
            plt.ylabel('Change in Free Energy (DelG)')
            plt.title('Temperature vs. DelG')
            plt.show()
            
        except:
            y4 = list(y4)
            xlst = []
            ylst = []
            for i in range(len(y4)):
                for j in range(len(y4[i])):
                    xlst.append(x4[i])
                    ylst.append(y4[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Temp (C)')
            plt.ylabel('DelG')
            plt.title('Temp vs. DelG Psat Curve')
            plt.legend(title = 'Psat Curve')
        plt.figure()
        # Plots for delV
        try: 
            pDelVPlot = sorted(pDelV.items()) # sorted by key, return a list of tuples
            x5, y5 = zip(*pDelVPlot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x5, y5)
            if self.psat == False:
                plt.legend(temper, title = "Temperatures (C)")
            plt.xlabel('Pressure (bar)')
            plt.ylabel('Change in Volume (DelV)')
            plt.title('Pressure vs. DelV')
            plt.show()
        except:
            y5 = list(y5)
            xlst =[]
            ylst = []
            for i in range(len(y5)):
                for j in range(len(y5[i])):
                    xlst.append(x5[i])
                    ylst.append(y5[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Temp (C)')
                    
            plt.xlabel('Pressure (bar)')
            plt.ylabel('DelV')
            plt.title('Pressure vs. DelV Psat Curve')
                    
        plt.figure()            
        try:
            tDelVPlot = sorted(tDelV.items()) # sorted by key, return a list of tuples
            x6, y6 = zip(*tDelVPlot) # unpack a list of pairs into two tuples
            plt.figure()
            plt.plot(x6, y6)
            plt.legend(press, title = "Pressure (Bar)")
            plt.xlabel('Temperature (C)')
            plt.ylabel('Change in Volume (DelV)')
            plt.title('Temperature vs. DelV')
            plt.show()
        except:
            xlst = []
            ylst = []
            y6 = list(y6)
            for i in range(len(y6)):
                for j in range(len(y6[i])):
                    xlst.append(x6[i])
                    ylst.append(y6[i][j])
            plt.plot(xlst,ylst)
            plt.xlabel('Temp (C)')
            plt.ylabel('DelV')
            plt.title('Temp vs. DelV Psat Curve')
        return
    
    def outLoop(self):
        '''A helper function to allow SUPCRTBL to run'''
        running = True
        while(running):
            line = self.pout.readline().decode(sys.stdout.encoding)
            print(line, end='')
            running='\n' in line
        print('Finished')
    
    def run_supcrt(self):
        '''A function that runs the pre-compiled SUPCRTBL found in the file folder'''
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
            filename = op.dirname(op.abspath(__file__)) + '\\resources\\' + customFile
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




