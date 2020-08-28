
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
import pkg_resources
import os.path as op

# ### Defining a Global Variables (Location and Constants)

# In[3]:
resource_package = 'DEWPython'
min_path = '/'.join(('resources', 'mineralDictionary.txt'))
mineralPath = pkg_resources.resource_filename(resource_package, min_path)
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
sup_path = '/'.join(('resources', 'SUPCRTBL.exe'))
supPath =  pkg_resources.resource_filename(resource_package, sup_path)

global Tr, bigQ, Chi, Pr, E_PrTr, bigR, Psi, Theta, Upsilon, Conversion, mineralDictionary
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
    

