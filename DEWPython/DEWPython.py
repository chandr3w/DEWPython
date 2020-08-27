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
from matplotlib.lines import Line2D
get_ipython().run_line_magic('matplotlib', 'inline')
from collections import defaultdict
import os.path as op

# ### Defining a Global Variables (Location and Constants)

# In[3]:
mineralPath = op.dirname(op.abspath(__file__)) + '/resources/mineralDictionary.txt'
gasPath = op.dirname(op.abspath(__file__)) + '/resources/gasLst.txt'
aqPath = op.dirname(op.abspath(__file__)) + '/resources/aqLst.txt'


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
        
        diaL = pd.read_csv('dielectric.csv', header = None)
        dia = diaL.to_numpy()
        dia = dia[4:, 1:]
        diaTrim = dia[1:, 1:]
        diaCollection = []
        for row in range(len(diaTrim)):
            for pressure in range(len(diaTrim[0])):
                # in form pressure, temperature, value
                diaCollection.append([dia[0][pressure + 1], dia[row + 1][0], diaTrim[row][pressure]])

        watDen = pd.read_csv('Wat_den.csv', header = None)
        w = watDen.to_numpy()
        w = w[4:, 1:]
        wTrim = w[1:,1:]
        watDenCollection = []
        for row in range(len(wTrim)):
            for pressure in range(len(wTrim[0])):
                # in form pressure, temperature, value
                watDenCollection.append([w[0][pressure + 1], w[row + 1][0], wTrim[row][pressure]])

        gibbsOfWater = pd.read_csv('water_gibbs.csv', header = None)
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
            ptSheet = pd.read_csv('input.csv',encoding= 'unicode_escape', header = None)
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
        self.proc = subprocess.Popen('SUPCRTBL.exe',shell = True, stdout = subprocess.PIPE, stdin = subprocess.PIPE, stderr = subprocess.STDOUT)
        self.pout = self.proc.stdout
        self.pin = self.proc.stdin
        threading.Thread(target=self.outLoop).start()
        while(self.proc.poll() is None):
            var = input('User Input: ')
            if '.txt' in var:
                self.supcrtFile = var
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
            filename = customFile
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


class DEWEquations:
    '''The class here imports all the equations that the authors of the Deep Earth Water Model Excel Sheet use 
    and converts them into Python'''
    def calculateDensity(pressure, temperature, equation, error, Psat):

        ''' Function to calculate the density of water. Essentially performs guesses and checks with
        different densities until it reaches the correct pressure down to two decimal places,
        as calculated by either Zhang & Duan (2005) or Zhang & Duan (2009).
        ---Input---
        pressure       - The pressure to calculate the density of water at, in bars
        temperature    - The temperature to calculate the density of water at, in Celsius
        equation       - Determines which equation of state to use in calculating the density.
                         equation = 1 corresponds to using Zhang & Duan (2005)
                         equation = 2 corresponds to using Zhang & Duan (2009)
        error          - This function uses a form of the bisection method. This variable indicates
                         how close the approximation should get. Eg. if error = 0.01, the density calculated
                         will calculate the pressure using the respective equation accurate to 0.01 of the input pressure
        Psat           - Determines if the polynomial fit to psat densities should be used in the event
                         that calculations are along the Psat curve
        ---Output---
        Returns the density of water at the input pressure and temperature, in units of g/cm^3. The density returned
        will calculate a pressure which differs from the input pressure by the value of "error" or less. If a proper value
        for the equation was not entered, zero is returned.
        '''
        fn_return_value = 0
        if Psat == True:

            #This equation models the density of water as a function of temperature along the Psat curve.
            #It has an R^2 value of 0.9999976885 as compared with Supcrt92 values.
            
            fn_return_value = - 1.01023381581205E-104 * pow(temperature, np.double(40)) + - 1.1368599785953E-27 * pow(temperature, np.double(10)) + - 2.11689207168779E-11 * pow(temperature, np.double(4)) + 1.26878850169523E-08 * pow(temperature, np.double(3)) + - 4.92010672693621E-06 * pow(temperature, np.double(2)) + - 3.2666598612692E-05 * temperature + 1.00046144613017
     
        else:
            #Define variables
            minGuess = 0.00001
            guess = 0.00001
            maxGuess = 7.5 * equation - 5
            calcP = 0
            #Loop through and find the density
            for i in range(1, 51):
                #Calculates the pressure using the specified equation
                calcP = DEWEquations.calculatePressure(guess, temperature, equation)
                #If the calculated pressure is not equal to input pressure, this determines a new
                #guess for the density based on current guess and how the calculated pressure
                #relates to the input pressure. In effect, this a form of a bisection method.
                if np.absolute(calcP - pressure) > error:
                    if calcP > pressure:
                        maxGuess = guess
                        guess = ( guess + minGuess )  / 2
                    elif calcP < pressure:
                        minGuess = guess
                        guess = ( guess + maxGuess )  / 2
                else:
                    fn_return_value = guess
                    break
        return fn_return_value
    
    

    def calculatePressure(density, temperature, equation):
        '''Calculates the pressure of water as a function of density and temperature using one of two
        equation of states.
        ---Input---
        density        - The density to use in finding a pressure, in g/cm^3
        temperature    - The temperature to use in finding a pressure, in Celsius
        equation       - The equation of state to use when calculating the pressure.
                         equation = 1 corresponds to using Zhang & Duan (2005)
                         equation = 2 corresponds to using Zhang & Duan (2009)
        ---Output---
        Returns the pressure of water corresponding to the input density and temperature, in units of bars.
        If a proper value for the equation was not entered, zero is returned.
        '''
        B = None

        C = None

        D = None

        E = None

        f = None

        g = None

        m = None
        m = np.double(18.01528)
        select_variable_0 = equation
        if (select_variable_0 == 1):
            ZD05_R = 83.144
            ZD05_Vc = 55.9480373
            ZD05_Tc = 647.25
            TK = temperature + 273.15
            Vr = m / density / ZD05_Vc
            Tr = TK / ZD05_Tc
            B = 0.349824207 - 2.91046273 /  ( Tr * Tr )  + 2.00914688 /  ( Tr * Tr * Tr )
            C = 0.112819964 + 0.748997714 /  ( Tr * Tr )  - 0.87320704 /  ( Tr * Tr * Tr )
            D = 0.0170609505 - 0.0146355822 /  ( Tr * Tr )  + 0.0579768283 /  ( Tr * Tr * Tr )
            E = - 0.000841246372 + 0.00495186474 /  ( Tr * Tr )  - 0.00916248538 /  ( Tr * Tr * Tr )
            f = - 0.100358152 / Tr
            g = np.double(- 0.00182674744 * Tr)
            delta = 1 + B / Vr + C /  ( Vr * Vr )  + D / pow(Vr, np.double(4)) + E / pow(Vr, np.double(5)) +  ( f /  ( Vr * Vr )  + g / pow(Vr, np.double(4)) )  * np.exp(- 0.0105999998 /  ( Vr * Vr ))
            fn_return_value = ZD05_R * TK * density * delta / m
        elif (select_variable_0 == 2):
            ZD09_R = 0.083145
            #Constant equal to ZD09_epsilon / (3.0626 * ZD09_omega^3)
            ZD09_c1 = 6.971118009
            #ZD09_epsilon = 510       'Lenard-Jones parameter in units of K
            #ZD09_omega = 2.88        'Lenard-Jones parameter in units of 1E-10 m
            #Prefactor calculated from 1000 * pow(ZD09_omega / 3.691, 3)
            dm = 475.05656886 * density
            #Prefactor calculated from 0.001 * pow(3.691 / ZD09_omega, 3)
            Vm = 0.0021050125 *  ( m / density )
            #Prefactor calculated from 154 / ZD09_epsilon
            Tm = 0.3019607843 *  ( temperature + 273.15 )  
            
            B = 0.029517729893 - 6337.56452413 /  ( Tm * Tm )  - 275265.428882 /  ( Tm * Tm * Tm )
            C = 0.00129128089283 - 145.797416153 /  ( Tm * Tm )  + 76593.8947237 /  ( Tm * Tm * Tm )
            D = 2.58661493537E-06 + 0.52126532146 /  ( Tm * Tm )  - 139.839523753 /  ( Tm * Tm * Tm )
            E = - 2.36335007175E-08 + 0.00535026383543 /  ( Tm * Tm )  - 0.27110649951 /  ( Tm * Tm * Tm )
            f = 25038.7836486 /  ( Tm * Tm * Tm )
            delta = 1 + B / Vm + C /  ( Vm * Vm )  + D / pow(Vm, 4) + E / pow(Vm, 5) + f /  ( Vm * Vm )  *  ( 0.73226726041 + 0.015483335997 /  ( Vm * Vm ) )  * np.exp(- 0.015483335997 /  ( Vm * Vm ))
            Pm = ZD09_R * Tm * delta / Vm
            fn_return_value = Pm * ZD09_c1
        else:
            fn_return_value = 0
        return fn_return_value

    
    
    
    def calculate_drhodP(density, temperature, equation):
        '''Calculates the partial derivative of density with respect to pressure, i.e. (d(rho)/dP)_T
        This is done using one of two equations of state for water.
        ---Input---
        density        - The density of water, in g/cm^3
        temperature    - The temperature of water, in Celsius
        equation       - The equation of state to use when calculating the pressure.
                         equation = 1 corresponds to using Zhang & Duan (2005)
                         equation = 2 corresponds to using Zhang & Duan (2009)
        ---Output---
        Returns the partial derivative of density with respect to pressure of water corresponding
        to the input density and temperature, in units of g^3/cm^3/bar. If a proper value for the equation
        was not entered, zero is returned.
        '''
        B = None

        C = None

        D = None

        E = None

        f = None

        g = None

        m = None
        m = np.double(18.01528)
        select_variable_1 = equation
        if (select_variable_1 == 1):
            ZD05_R = 83.144
            ZD05_Vc = 55.9480373
            ZD05_Tc = 647.25
            TK = np.double(temperature + 273.15)
            Tr = TK / ZD05_Tc
            cc = ZD05_Vc / m
            Vr = m /  ( density * ZD05_Vc )
            B = 0.349824207 - 2.91046273 /  ( Tr * Tr )  + 2.00914688 /  ( Tr * Tr * Tr )
            C = 0.112819964 + 0.748997714 /  ( Tr * Tr )  - 0.87320704 /  ( Tr * Tr * Tr )
            D = 0.0170609505 - 0.0146355822 /  ( Tr * Tr )  + 0.0579768283 /  ( Tr * Tr * Tr )
            E = - 0.000841246372 + 0.00495186474 /  ( Tr * Tr )  - 0.00916248538 /  ( Tr * Tr * Tr )
            f = - 0.100358152 / Tr
            g = np.double(- 0.00182674744 * Tr)
            delta = 1 + B / Vr + C /  ( Vr * Vr )  + D / pow(Vr, 4) + E / pow(Vr, 5) +  ( f /  ( Vr * Vr )  + g / pow(Vr, 4) )  * np.exp(- 0.0105999998 / pow(Vr, 2))
            kappa = B * cc + 2 * C *  ( cc * cc )  * density + 4 * D * pow(cc, 4) * pow(density, 3) + 5 * E * pow(cc, 5) * pow(density, 4) +  ( 2 * f *  ( cc * cc )  * density + 4 * g * pow(cc, 4) * pow(density, 3) -  ( f /  ( Vr * Vr )  + g / pow(Vr, 4) )  *  ( 2 * 0.0105999998 *  ( cc * cc )  * density ) )  * np.exp(- 0.0105999998 /  ( Vr * Vr ))
            fn_return_value = m /  ( ZD05_R * TK *  ( delta + density * kappa ) )
        elif (select_variable_1 == 2):
            ZD09_R = 0.083145
            ZD09_c1 = 6.971118009
            #ZD09_epsilon = 510       'Lenard-Jones parameter in units of K
            #ZD09_omega = 2.88        'Lenard-Jones parameter in units of 1E-10 m
            #Prefactor calculated from 1000 * pow(ZD09_omega / 3.691, 3)
            dm = 475.05656886 * density
            #Prefactor calculated from 0.001 * pow(3.691 / ZD09_omega, 3)
            Vm = 0.0021050125 *  ( m / density )
            #Prefactor calculated from 154 / ZD09_epsilon
            Tm = 0.3019607843 *  ( temperature + 273.15 )   
            B = 0.029517729893 - 6337.56452413 /  ( Tm * Tm )  - 275265.428882 /  ( Tm * Tm * Tm )
            C = 0.00129128089283 - 145.797416153 /  ( Tm * Tm )  + 76593.8947237 /  ( Tm * Tm * Tm )
            D = 2.58661493537E-06 + 0.52126532146 /  ( Tm * Tm )  - 139.839523753 /  ( Tm * Tm * Tm )
            E = - 2.36335007175E-08 + 0.00535026383543 /  ( Tm * Tm )  - 0.27110649951 /  ( Tm * Tm * Tm )
            f = 25038.7836486 /  ( Tm * Tm * Tm )
            delta = 1 + B / Vm + C /  ( Vm * Vm )  + D / pow(Vm, 4) + E / pow(Vm, 5) + f /  ( Vm * Vm )  *  ( 0.73226726041 + 0.015483335997 /  ( Vm * Vm ) )  * np.exp(- 0.015483335997 /  ( Vm * Vm ))
            kappa = B / m + 2 * C * dm /  ( m * m )  + 4 * D * pow(dm, 3) / pow(m, 4) + 5 * E * pow(dm, 4) / pow(m, 5) + ( 2 * f * dm /  ( m * m )  *  ( 0.73226726041 + 0.015483335997 /  ( Vm * Vm ) )  + f / pow(Vm, 2) *  ( 1 - 0.73226726041 - 0.015483335997 /  ( Vm * Vm ) )  *  ( 2 * 0.015483335997 * dm /  ( m * m ) ) )  * np.exp(- 0.015483335997 /  ( Vm * Vm ))
            
            ##### Adding  a comment here because I've made ZD09_c4 into ZD09 C_1 #######
            ##### Original line######
            #fn_return_value = ZD09_c1 * m /  ( ZD09_c4 * ZD09_R * Tm *  ( delta + dm * kappa ) )
            fn_return_value = ZD09_c1 * m /  ( ZD09_c1 * ZD09_R * Tm *  ( delta + dm * kappa ) )
        else:
            fn_return_value = 0
        return fn_return_value
    
    
    
    
    def calculateGibbsOfWater(pressure, temp, equation, densityEquation, Psat):
        '''This function calculates the Gibbs Free Energy of Water. It can calculate with two equations.
        ---Input---'
        pressure           - The pressure to calculate the Gibbs Free Energy at, in bars
        temperature        - The temperature to calculate the Gibbs Free Energy at, in Celsius
        equation           - Determines which equation to use to calculate the Gibbs Free Energy,
                             either Delaney & Helgeson (1978), corresonding to equation = 1, or simply integrating
                             over the volume of water, corresponding to equation = 2
        density Equation    - Determines which equation to use to find the density, and thus the volume of water.
        Psat               - Determines if the calculation should be done at Psat.
        ---Output---
        Returns the Gibbs Free Energy of water in units of cal/mol. If a proper value for equation was not entered,
        zero is returned.
        '''
        if Psat == True:
            #This equation models the Gibbs Free Energy of water as a function of temperature along the Psat curve.
            #It has an R^2 value of 0.9999999984518 as compared with Supcrt92 values.
            fn_return_value = - 2.72980941772081E-103 * pow(temp, np.double(40)) + 2.88918186300446E-25 * pow(temp, np.double(10)) + - 2.21891314234246E-08 * pow(temp, np.double(4)) + 3.0912103873633E-05 * pow(temp, np.double(3)) + - 3.20873264480928E-02 * pow(temp, np.double(2)) + - 15.169458452209 * temp + - 56289.0379433809
        else:
            select_variable_2 = equation
            if (select_variable_2 == 1):
                coeff = {}
                coeff[0] = - 56130.073
                coeff[1] = 0.38101798
                coeff[2] = - 0.0000021167697
                coeff[3] = 2.0266445E-11
                coeff[4] = - 8.3225572E-17
                coeff[5] = - 15.285559
                coeff[6] = 0.0001375239
                coeff[7] = - 1.5586868E-09
                coeff[8] = 6.6329577E-15
                coeff[9] = - 0.026092451
                coeff[10] = 0.000000035988857
                coeff[11] = - 2.7916588E-14
                coeff[12] = 0.000017140501
                coeff[13] = - 1.6860893E-11
                coeff[14] = - 6.0126987E-09
                gibbsFreeEnergy = 0
                Count = 0
                
                for j in range(0, 5):
                    for k in range(0, 5 - j):
                        temp = np.absolute(temp)

                        gibbsFreeEnergy = gibbsFreeEnergy + coeff[Count] * pow((temp), np.double(j)) * pow(pressure, np.double(k))
                        
                        Count = Count + 1
                fn_return_value = gibbsFreeEnergy
            elif (select_variable_2 == 2):
                
                #then defines the gibbs free energy as the integral over the volume as a function of temperature.
                #We can only perform this calculation if we can use one of the two density equations included
                #in the code. If densityEquation equals three, then that implies the user chose to use custom
                #density values. Because this procedure requires integration over a range of densities, this
                #cannot be calculated if the user has custom density values. Therefore, this will just return zero.
                if ( densityEquation == 3 ) :
                    fn_return_value = 0
                    
                #Gibbs Free Energy of water at 1 kb. This equation is a polynomial fit to data as a function of temperature.
                #It is valid in the range of 100 to 1000 C.

                temp = np.absolute(temp) 
                GAtOneKb = 2.6880734E-09 *(temp * temp)*(temp*temp) + 0.00000063163061 * (temp * temp * temp) - 0.019372355 *  ( temp * temp )  - 16.945093 * temp - 55769.287
                
                
                if pressure < 1000:
                    fn_return_value = 0
                elif pressure == 1000:
                    fn_return_value = GAtOneKb
                elif pressure > 1000:
                    integral = 0
                    #Integral is sum of rectangles with this width. This function in effect limits the spacing
                    #to 20 bars so that very small pressures do not have unreasonably small widths. Otherwise the width
                    #is chosen such that there are always 500 steps in the numerical integration. This ensures that for very
                    #high pressures, there are not a huge number of steps calculated which is very computationally taxing.
                    if ( pressure - 1000 )  / 500 < 20:
                        spacing = 20
                    else: 
                        spacing = ( pressure - 1000 )  / 500
                    
                    for i in range(1000, pressure + 1, spacing):
                        #This integral determines the density only down to an error of 100 bars
                        #rather than the standard of 0.01. This is done to save computational
                        #time. Tests indicate this reduces the computation by about a half while
                        #introducing little error from the standard of 0.01.
                        
                        integral = integral +  ( 18.01528 / DEWEquations.calculateDensity(i, temp, densityEquation, 100, False) / 41.84 )  * spacing
                        
                    fn_return_value = GAtOneKb + integral
                    
            else:
                fn_return_value = 0
        return fn_return_value
    
    
    
    
    def calculateEpsilon(density, temperature, equation, Psat):
        ''' This function calculates the dielectric constant (epsilon) of water using one of four possible equations.
        ---Input---
        density        - The density of water to use in calculating epsilon, in g/cm^3
        temperature    - The temperature to calculate epsilon with, in Celsius
        equation       - Determines which equation should be used to calculate the dielectric constant of water.
                         equation = 1 corresponds to using Johnson & Norton (1991), the equation used in Supcrt
                         equation = 2 corresponds to using Franck (1990)
                         equation = 3 corresponds to using Fernandez (1997)
                         equation = 4 corredponds to using the Power Function. This is an equation derived by
                         Dimitri Sverjensky and Brandon Harison at Johns Hopkins University.
        Psat           - Determines if the polynomial fit to psat dielectric constant values should be used
                         in the event that calculations are along the Psat curve
        ---Output---
        Returns the Dielectric constant of water at the given density and temperature. If a proper value
        for equation was not entered, zero is returned.
        '''
        if Psat == True:
            #This equation models the dielectric constant of water as a function of temperature along the Psat curve.
            #It has an R^2 value of 0.9999991719 as compared with Supcrt92 values.
            fn_return_value = - 1.66686763214295E-77 * pow(temperature, np.double(30)) + - 9.02887020379887E-07 * pow(temperature, np.double(3)) + 8.4590281449009E-04 * pow(temperature, np.double(2)) + - 0.396542037778945 * temperature + 87.605024245432
        else:
            select_variable_3 = equation
            if (select_variable_3 == 1):
                T_hat = ( temperature + 273.15 )  / 298.15
                k0 = 1
                k1 = 14.70333593 / T_hat
                k2 = 212.8462733 / T_hat - 115.4445173 + 19.55210915 * T_hat
                k3 = - 83.3034798 / T_hat + 32.13240048 * T_hat - 6.69409865 *  ( T_hat * T_hat )
                k4 = - 37.86202045 /  ( T_hat * T_hat )  + 68.87359646 / T_hat - 27.29401652
                fn_return_value = k0 + k1 * density + k2 *  ( density * density )  + k3 * pow(density, 3) + k4 * pow(density, 4)
            elif (select_variable_3 == 2):
                pi = 3.14159265358979
                omega = 0.0000000268
                k = 1.380648E-16
                Na = 6.022E+23
                mu = 2.33E-18
                rhostar = ( density * 0.055508 )  * pow(omega, 3) * Na
                mustarsq = pow(mu, 2) /  ( k *  ( temperature + 273.15 )  * pow(omega, 3) )
                y = ( 4 * pi / 9 )  * rhostar * mustarsq
                f1 = 0.4341 * pow(rhostar, 2)
                f2 = - ( 0.05 + 0.75 * pow(rhostar, 3) )
                f3 = - 0.026 * pow(rhostar, 2) + 0.173 * pow(rhostar, 4)
                fn_return_value = ( ( 3 * y )  /  ( 1 - f1 * y ) )  *  ( 1 +  ( 1 - f1 )  * y + f2 *  ( y * y )  + f3 *  ( y * y * y ) )  + 1
            elif (select_variable_3 == 3):
                #Values for N_k
                N_k = {}
                N_k[0] = 0.978224486826
                N_k[1] = - 0.957771379375
                N_k[2] = 0.237511794148
                N_k[3] = 0.714692224396
                N_k[4] = - 0.298217036956
                N_k[5] = - 0.108863472196
                N_k[6] = 0.0949327488264
                N_k[7] = - 0.00980469816509
                N_k[8] = 0.000016516763497
                N_k[9] = 9.37359795772E-05
                N_k[10] = - 1.2317921872E-10
                N_k[11] = 0.00196096504426
                #Values for i_k
                i_k = {}
                i_k[0] = 1
                i_k[1] = 1
                i_k[2] = 1
                i_k[3] = 2
                i_k[4] = 3
                i_k[5] = 3
                i_k[6] = 4
                i_k[7] = 5
                i_k[8] = 6
                i_k[9] = 7
                i_k[10] = 10
                #Values for j_k
                j_k = {}
                j_k[0] = 0.25
                j_k[1] = 1
                j_k[2] = 2.5
                j_k[3] = 1.5
                j_k[4] = 1.5
                j_k[5] = 2.5
                j_k[6] = 2
                j_k[7] = 2
                j_k[8] = 5
                j_k[9] = 0.5
                j_k[10] = 10
                avogadro = 6.0221367E+23
                dipole = 6.138E-30
                epsilon_o = 8.8541878176204E-12
                boltzmann = 1.380658E-23
                alpha = 1.636E-40
                density_c = 17873.728
                T_c = 647.096
                #Convert density and temperature units
                density_molm3 = density * 0.055508 * 1000000
                T_K = temperature + 273.15
                #Defining the g equation
                g = 1
                for ii in range(0, 11):
                    g = g + N_k[ii] * pow(density_molm3 / density_c, np.double(i_k[ii])) * pow(T_c / T_K, np.double(j_k[ii]))
                g = g + N_k[11] *  ( density_molm3 / density_c )  * pow(T_K / 228 - 1, - 1.2)
                #Defining the A, B, and C equations
                A = ( avogadro * pow(dipole, 2) * density_molm3 * g )  /  ( epsilon_o * boltzmann * T_K )
                B = ( avogadro * alpha * density_molm3 )  /  ( 3 * epsilon_o )
                C = 9 + 2 * A + 18 * B + A * A + 10 * A * B + 9 * B * B
                fn_return_value = ( 1 + A + 5 * B + np.sqrt(C) )  /  ( 4 - 4 * B )
            elif (select_variable_3 == 4):
                #Relevant parameters
                a1 = - 1.57637700752506E-03
                a2 = 6.81028783422197E-02
                a3 = 0.754875480393944
                b1 = - 8.01665106535394E-05
                b2 = - 6.87161761831994E-02
                b3 = 4.74797272182151
                A = a1 * temperature + a2 * np.sqrt(temperature) + a3
                B = b1 * temperature + b2 * np.sqrt(temperature) + b3
                fn_return_value = np.exp(B) * pow(density, np.double(A))
            else:
                fn_return_value = 0
        return fn_return_value
    
    
    
    
    def calculate_depsdrho(density, temperature, equation):
        '''Calculates the partial derivative of the dielectric constant (epsilon) with respect to density, i.e. (d(eps)/d(rho))_T
        This is done using one of four possible equations
        ---Input---
        density        - The density of water to calculate with, in g/cm^3
        temperature    - The temperature to calculate with, in Celsius
        equation       - Determines which equation should be used to calculate the derivative
                         equation = 1 corresponds to using Johnson & Norton (1991), the equation used in Supcrt
                         equation = 2 corresponds to using Franck (1990)
                         equation = 3 corresponds to using Fernandez (1997)
                         equation = 4 corredponds to using the Power Function. This is an equation derived by
                         Dimitri Sverjensky and Brandon Harison at Johns Hopkins University.
        ---Output---
        Returns the partial derivative of the dielectric constant with respect to density in units of cm^3/g. If a proper value
        for equation was not entered, zero is returned.
        '''
        select_variable_4 = equation
        if (select_variable_4 == 1):
            T_hat = ( temperature + 273.15 )  / 298.15
            k1 = 14.70333593 / T_hat
            k2 = 212.8462733 / T_hat - 115.4445173 + 19.55210915 * T_hat
            k3 = - 83.3034798 / T_hat + 32.13240048 * T_hat - 6.69409865 *  ( T_hat * T_hat )
            k4 = - 37.86202045 /  ( T_hat * T_hat )  + 68.87359646 / T_hat - 27.29401652
            fn_return_value = k1 + 2 * k2 * density + 3 * k3 * pow(density, 2) + 4 * k4 * pow(density, 3)
        elif (select_variable_4 == 2):
            pi = 3.14159265358979
            omega = 0.0000000268
            k = 1.380648E-16
            Na = 6.022E+23
            mu = 2.33E-18
            density = density * 0.055508
            cc = pow(omega, 3) * Na
            rhostar = density * cc
            mustarsq = pow(mu, 2) /  ( k *  ( temperature + 273.15 )  * pow(omega, 3) )
            y = ( 4 * pi / 9 )  * rhostar * mustarsq
            f1 = 0.4341 * pow(rhostar, 2)
            f2 = - ( 0.05 + 0.75 * pow(rhostar, 3) )
            f3 = - 0.026 * pow(rhostar, 2) + 0.173 * pow(rhostar, 4)
            dydrho = ( 4 * pi / 9 )  * mustarsq * cc
            df1drho = 2 * 0.4341 * pow(cc, 2) * density
            df2drho = - 3 * 0.75 * pow(cc, 3) * pow(density, 2)
            df3drho = - 2 * 0.026 * pow(cc, 2) * density + 4 * 0.173 * pow(cc, 4) * pow(density, 3)
            eps = ( ( 3 * y )  /  ( 1 - f1 * y ) )  *  ( 1 +  ( 1 - f1 )  * y + f2 *  ( y * y )  + f3 *  ( y * y * y ) )  + 1
            #The 0.055508 value converts the units from cm^3/mol to cm^3/g
            fn_return_value = 0.05508 *  ( ( ( dydrho + pow(y, 2) * df1drho )  /  ( 1 - f1 * y ) )  *  ( eps - 1 )  / y +  ( ( 3 * y )  /  ( 1 - f1 * y ) )  *  
                                          ( - df1drho * y + df2drho * pow(y, 2) + df3drho * pow(y, 3) +  ( 1 - f1 + 2 * f2 * y + 3 * f3 * y * y )  * dydrho ) )
        elif (select_variable_4 == 3):
            #Values for N_k
            N_k = {}
            N_k[0] = 0.978224486826
            N_k[1] = - 0.957771379375
            N_k[2] = 0.237511794148
            N_k[3] = 0.714692224396
            N_k[4] = - 0.298217036956
            N_k[5] = - 0.108863472196
            N_k[6] = 0.0949327488264
            N_k[7] = - 0.00980469816509
            N_k[8] = 0.000016516763497
            N_k[9] = 9.37359795772E-05
            N_k[10] = - 1.2317921872E-10
            N_k[11] = 0.00196096504426
            #Values for i_k
            i_k = {}
            i_k[0] = 1
            i_k[1] = 1
            i_k[2] = 1
            i_k[3] = 2
            i_k[4] = 3
            i_k[5] = 3
            i_k[6] = 4
            i_k[7] = 5
            i_k[8] = 6
            i_k[9] = 7
            i_k[10] = 10
            #Values for j_k
            j_k = {}
            j_k[0] = 0.25
            j_k[1] = 1
            j_k[2] = 2.5
            j_k[3] = 1.5
            j_k[4] = 1.5
            j_k[5] = 2.5
            j_k[6] = 2
            j_k[7] = 2
            j_k[8] = 5
            j_k[9] = 0.5
            j_k[10] = 10
            avogadro = 6.0221367E+23
            dipole = 6.138E-30
            epsilon_o = 8.8541878176204E-12
            boltzmann = 1.380658E-23
            alpha = 1.636E-40
            density_c = 17873.728
            T_c = 647.096
            #Convert density and temperature units
            density_molm3 = density * 0.055508 * 1000000
            T_K = temperature + 273.15
            #Defining the g equation
            g = 1
            for ii in range(0, 11):
                g = g + N_k[ii] * pow(density_molm3 / density_c, np.double(i_k[ii])) * pow(T_c / T_K, np.double(j_k[ii]))
            g = g + N_k[11] *  ( density_molm3 / density_c )  * pow(T_K / 228 - 1, - 1.2)
            #Defining the dgdrho equation
            dgdrho = 0
            for ii in range(0, 11):
                dgdrho = dgdrho + i_k[ii] * N_k[ii] *  ( pow(density_molm3, np.double(i_k[ii] - 1)) / pow(density_c, np.double(i_k[ii])) )  * pow(T_c / T_K, np.double(j_k[ii]))
            dgdrho = dgdrho +  ( N_k[11] / density_c )  * pow(T_K / 228 - 1, - 1.2)
            #Defining the A, B, and C equations
            A = ( avogadro * pow(dipole, 2) * density_molm3 * g )  /  ( epsilon_o * boltzmann * T_K )
            B = ( avogadro * alpha * density_molm3 )  /  ( 3 * epsilon_o )
            C = 9 + 2 * A + 18 * B + A * A + 10 * A * B + 9 * B * B
            #Defining the derivatives and epsilon
            dAdrho = A / density_molm3 +  ( A / g )  * dgdrho
            dBdrho = B / density_molm3
            dCdrho = 2 * dAdrho + 18 * dBdrho + 2 * A * dAdrho + 10 *  ( dAdrho * B + A * dBdrho )  + 18 * B * dBdrho
            eps = ( 1 + A + 5 * B + pow(np.double(C), 0.5))   /  ( 4 - 4 * B )
            #The 55508 value converts the units from m^3/mol to cm^3/g
            fn_return_value = 55508 *  ( 1 /  ( 4 - 4 * B ) )  *  ( 4 * dBdrho * eps + dAdrho + 5 * dBdrho + 0.5 * pow(np.double(C), - 0.5) * dCdrho )
        elif (select_variable_4 == 4):
            #Relevant parameters
            a1 = - 1.57637700752506E-03
            a2 = 6.81028783422197E-02
            a3 = 0.754875480393944
            b1 = - 8.01665106535394E-05
            b2 = - 6.87161761831994E-02
            b3 = 4.74797272182151
            A = a1 * temperature + a2 * np.sqrt(temperature) + a3
            B = b1 * temperature + b2 * np.sqrt(temperature) + b3
            fn_return_value = A * np.exp(B) * pow(density, A - 1)
        else:
            fn_return_value = 0
        return fn_return_value
    
    
    
    def calculateOmega(P, T, density, name, wref, Z):
        '''This function calculates the born coefficient omega for aqueous species as a function of pressure and temeprature
        ---Input---
        P          - Pressure to calculate at, in bars
        T          - Temperature to calculate at, in Celsius
        density    - Density of water to calculate omega at, in g/cm^3. This could be calculated from P and T, but
                     it is used as an input parameter to save on calculation time.
        name       - The name of the species this is being calculated for.
        wref       - The value of omega at standard pressure and temperature, in units of cal/mol. This should not be
                     the value generally given as omega*1E-5, but rather the actual value of omega.
        Z          - The charge of the species
        ---Output---
        Returns the value of omega at the input P and T. If Z is zero, the wprtr value is used. The value returned is
        in units of cal/mol and NOT multiplied by 10^-5.
        '''
        #If species is hydrogen, the species is neutral, or the pressure is above 6 kb,
        #this equation is not necessary because omega is very close to wref.
        if name == 'H+' or Z == 0 or P > 6000:
            fn_return_value = wref
        else:
            #These equations are given by Shock et al. (1992)
            eta = 166027
            #Defines the electrostatic radius at reference pressure and temperature
            reref = Z * Z /  ( wref / eta + Z / 3.082 )
            #This represents the pressure and temperature dependent solvent function
            g = DEWEquations.calculateG(P, T, density)
            #Defines the electrostatic radius at the input P and T
            re = reref + (Z) * g
            fn_return_value = eta *  ( Z * Z / re - Z /  ( 3.082 + g ) )
        return fn_return_value
    
    
    
    def calculateG(P, T, density):
        '''Calculates the pressure and temperature dependent solvent function. This function should only be
        used for pressures less than 6 kb.
        ---Input---
        P          - The pressure to calculate at, in bars
        T          - The temperature to calculate at, in celsius
        density    - The density of water at which to calculate g at, in g/cm^3
        ---Output---
        Returns the value of the g function. If the density is greather than 1 g/cm^3, then zero is returned.'''
        if density >= 1:
            fn_return_value = 0
        else:
            a_g = - 2.037662 + 0.005747 * T - 0.000006557892 * T * T
            b_g = 6.107361 - 0.01074377 * T + 0.00001268348 * T * T
            #Calculates the difference function in the case where we need to calculate at Psat conditions
            if ( P <= 1000 and T >= 155 and T <= 355 ) :
                f = ( pow(( T - 155 )  / 300, 4.8) + 36.66666 * pow(( T - 155 )  / 300, np.double(16)) )  *( - 1.504956E-10 * pow(1000 - P, np.double(3)) + 5.017997E-14 * pow(1000 - P, np.double(4)) )
            else:
                f = 0
            fn_return_value = a_g * pow(1 - density, b_g) - f
        return fn_return_value
    
    def calculate_domegadP(P, T, density, name, wref, Z, densityEquation, Psat):
        '''This function calculates the derivative of the born coefficient omega with respect to pressure
        for aqueous species as a function of pressure and temeprature
        ---Input---
        P                  - Pressure to calculate at, in bars
        T                  - Temperature to calculate at, in Celsius
        density            - Density of water to calculate omega at, in g/cm^3. This could be calculated from P and T, but
                             it is used as an input parameter to save on calculation time.
        name               - The name of the species this is being calculated for.
        wref               - The value of omega at standard pressure and temperature, in units of cal/mol. This should not be
                             the value generally given as omega*1E-5, but rather the actual value of omega.
        Z                  - The charge of the species
        densityEquation    - Determines which equation to use in calculating the derivative of density
                             with respect to pressure. This is passed direction to calculate_dgdP
                             equation = 1  corresponds to Zhang & Duan (2005)
                             equation = 1  corresponds to Zhang & Duan (2009)
        Psat               - Determines if the calculation should be done along the Psat curve. In this case
                             there is no equation for drhodP and a polynomial fit to data from Shock et al. (1992) is used.
        ---Output---
        Returns the value of the derivative of omega with respect to pressure at the input P and T. If Z is zero, then
        the derivative is zero. The value returned is in units of cal/mol/bar
        '''
        #If species is hydrogen, the species is neutral, or the pressure is above 6 kb,
        #this equation is not necessary because omega is very close to wref.
        if name == 'H+' or Z == 0 or P > 6000:
            fn_return_value = 0
        else:
            #These equations are given by Shock et al. (1992)
            eta = 166027
            #Defines the electrostatic radius at reference pressure and temperature
            reref = Z * Z /  ( wref / eta + Z / 3.082 )
            #This represents the pressure and temperature dependent solvent function and its derivative
            g = DEWEquations.calculateG(P, T, density)
            dgdP = DEWEquations.calculate_dgdP(P, T, density, g, densityEquation, Psat)
            #Defines the electrostatic radius at the input P and T
            re = reref + np.absolute(Z) * g
            fn_return_value = - eta *  ( np.absolute(Z * Z * Z) / pow(re, 2) - Z / pow(3.082 + g, 2) )  * dgdP
        return fn_return_value
    
    
    def calculate_dgdP(P, T, density, g, equation, Psat = True):
        '''Calculates the pressure derivative of the pressure and temperature dependent solvent function.
        This function should only be used for pressures less than 6 kb.
        ---Input---
        P          - The pressure to calculate at, in bars
        T          - The temperature to calculate at, in celsius
        density    - The density of water at which to calculate g at, in g/cm^3
        g          - The value of the g solvent function at the input P and T
        equation   - Determines which equation to use in calculating the derivative of density
                     with respect to pressure
                     equation = 1  corresponds to Zhang & Duan (2005)
                     equation = 1  corresponds to Zhang & Duan (2009)
        Psat       - Determines if the calculation should be done along the Psat curve. In this case
                     there is no equation for drhodP and a polynomial fit to data from Shock et al. (1992) is used.
        ---Output---
        Returns the pressure derivative of the g function. If the density is greather than 1 g/cm^3, then zero is returned.
        '''
        if Psat == True:
            #This equation models the derivative of the g solvent function with respect to pressure and
            #as a function of temperature along the Psat curve.
            #It has an R^2 value of 0.99995027718 as compared with values listed in Shock et al. (1992).
            #Particular care was taken to properly model the values at low temperatures which is why this
            #function not simply a polynomial
            if T < 0.01:
                fn_return_value = 0
            else:
                fn_return_value = np.exp(1.37105493109451E-10 * pow(np.log(T), np.double(15)) + - 1.43605469318795E-06 * pow(np.log(T), np.double(10)) + 26.2649453651117 * np.log(T) + - 125.108856715714) * 0.000001
        else:
            if density >= 1:
                fn_return_value = 0
            else:
                b_g = 6.107361 - 0.01074377 * T + 0.00001268348 * T * T
                #Calculates the difference function in the case where we need to calculate at Psat conditions
                if ( P <= 1000 and T >= 155 and T <= 355 ) :
                    dfdP = - ( pow(( T - 155 )  / 300, 4.8) + 36.66666 * pow(( T - 155 )  / 300, 16) )  *  ( 3 * - 1.504956E-10 * pow(1000 - P, 2) + 4 * 5.017997E-14 * pow(1000 - P, 3) )
                else:
                    dfdP = 0
                fn_return_value = - b_g * calculate_drhodP(density, T, equation) * g /  ( 1 - density )  - dfdP
        return fn_return_value
    
    def calculateQ(pressure, temperature, density, densityEquation, epsilonEquation, Psat):
        '''This method calculates the Born Coefficient Q as (1/eps^2)*(d(eps)/dP) - In other words the derivative of
        epsilon with respect to pressure, divided by epsilon squared
        ---Input---
        pressure           - The pressure to calculate Q at, in bars
        temperature        - The temperature to calculate Q at, in Celsius
        density            - The density at the input pressure and temperature, input simply to save time, in g/cm^3
        denistyEquation    - The density equation to use in calculating the density of water.
        epsilonEquation    - The epsilon equation to use in calculating epsilon.
        Psat               - Determines if the calculation should be done at Psat.
        ---Output---
        Outputs the value of Q in units of bar^-1
        Calculates the pressure and temperature dependent solvent function. This function should only be
        used for pressures less than 6 kb.
        ---Input---
        P          - The pressure to calculate at, in bars
        T          - The temperature to calculate at, in celsius
        density    - The density of water at which to calculate g at, in g/cm^3
        ---Output---
        Returns the value of the g function. If the density is greather than 1 g/cm^3, then zero is returned.
        '''
        if Psat == True:
            #This equation models the Q Born Coefficent as a function of temperature along the Psat curve.
            #It has an R^2 value of 0.99999998602 as compared with values listed in Shock et al. (1992).
            fn_return_value = ( 1.99258688758345E-49 * pow(temperature, np.double(20)) + - 4.43690270750774E-14 * pow(temperature, np.double(6)) + 4.29110215680165E-11 * pow(temperature, np.double(5)) + - 1.07146606081182E-08 * pow(temperature, np.double(4)) + 1.09982931856694E-06 * pow(temperature, np.double(3)) + 9.60705240954956E-06 * pow(temperature, np.double(2)) + 0.642579832259358 )  * 0.000001
        else:
            #This commented section is the code to calculate the value of Q using a finite difference derivative.
            #-------------------------
            #        Dim epsilon, delta, epsilonPlusDelta As Double
            #
            #        delta = 1
            #
            #        epsilon = DEWEquations.calculateEpsilon(density, temperature, epsilonEquation, False)
            #
            #        epsilonPlusDelta = DEWEquations.calculateEpsilon(calculateDensity(pressure + delta, temperature, densityEquation, 0.01, False), temperature, epsilonEquation, False)
            #
            #        calculateQ = (1 / pow(np.double(epsilon), 2)) * ((epsilonPlusDelta - epsilon) / delta)
            #-------------------------
            eps = DEWEquations.calculateEpsilon(density, temperature, epsilonEquation, Psat)
            depsdrho = DEWEquations.calculate_depsdrho(density, temperature, epsilonEquation)
            drhodP = DEWEquations.calculate_drhodP(density, temperature, densityEquation)
            fn_return_value = depsdrho * drhodP /  ( eps * eps )
        return fn_return_value
    


# In[ ]:




