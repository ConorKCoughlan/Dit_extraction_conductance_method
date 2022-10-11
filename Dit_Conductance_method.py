#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import glob
import matplotlib.pyplot as plt
import numpy as np
import re
import os


#Define frequencies (in terms of omega = 2*pi*f) that conductance was measured at for each temperature
freqvalues_omega = 2*np.pi*np.array([1e3,2e3,4e3,5e3,6e3,8e3,10e3,12e3,14e3,16e3,20e3,25e3,30e3,40e3,50e3,75e3,
                                      100e3,150e3,200e3,300e3,400e3,500e3,600e3,800e3,1e6]) 


#define voltage index (i.e. voltage datapoints in GV/CV measurements) in conductance peak of G/omega vs omega is shown, this involves a bit of trial and error until you can see clear peaks in the middle of G/omega vs omega plot
#e.g. Voltage indices at 200K, 125K and 77K
Voltage_index_200 = range(286,294,1) #200K
Voltage_index_125 = range(254,288,2) #125K
Voltage_index_77 = range(160,190,2) #77K

#constants
k_B = 1.38e-23 #Boltzmann's constant
m_0 = 9.1e-31 #electron mass
h = 6.626e-34 #Planck's constant
sigma = 1e-15 #trap capture cross-section, for Ge this can vary from 1e-13 to 1e-17. Worth looking up approximate literature values
q = 1.6e-19 #electron charge
epsilon_0 = 8.85e-12 #permittivity of free space
epsilon_Ge = 16.2 #Ge relative permittivity - can adjust depending on what semiconductor material is in MOSCAP

A = np.pi*(500e-6)**2/4 #Device area (m^2), currently assuming 500um size

class Ditextract:
    def __init__(self,directory,filenames_extracted,temperature,omegavalues,voltage_array):
        'Initialises variables within individual class'
        self.directory = directory #Directory where raw data files are located
        self.filenames_extracted = filenames_extracted #What filenames within directory are extracted
        self.temperature = temperature #Measurement temperature
        self.omegavalues = omegavalues #Refers to 'freqvalues_omega' above
        self.voltage_array = voltage_array #Refers to 'Voltage_index' above
        
        
    def gotodirectory(self):
        'Moves to relevant directory'
        os.chdir(self.directory)
        
        
    def extractfilenames(self,filenames_extracted):
        'Extracts files (typically .csv files) from directory'
        filenames = glob.glob(filenames_extracted)

        def atoi(text):
            'Converts digits in .csv files from strings to int - handy for listing files in order of measured frequency'
            return int(text) if text.isdigit() else text

        def natural_keys(text):
            'key used in .sort function below for sorting in order of measured frequency'
            return [ atoi(c) for c in re.split(r'(\d+)', text) ]

        filenames.sort(key  = natural_keys)
        
        return filenames
        
    def extract_parameters(self, filename):
        'Returns experimental parameters from .csv files - these are Voltage, Conductance, Cqpacitance and Frequency'
        listnum = [[]]  
        with open(filename) as csvfile:
            lines = csvfile.readlines()
            
            for line in lines:
                if 'DataValue' in line:
                    valuefloat = []
                    values = line.split(',')[1:]
                    for value in values:
                        valuefloat.append(float(value))
                        
                    listnum.append(valuefloat)
                    
                    
            del listnum[0]
                    
        
        
        Voltage = []
        Conductance = []
        Capacitance = []
        InvC = []
        
        
        for i in range(len(listnum)):
            
            Voltage.append(listnum[i][0])
            Conductance.append(listnum[i][3])
            Capacitance.append(listnum[i][2])
            InvC.append(1/(listnum[i][2])**2) #This outputs 1/C^2, handy for determining doping density of semiconductor
        
        return Voltage, Conductance, Capacitance, InvC
    
    
    
    def extract_parameters_alt(self,filename):
        'Returns same experimental parameters, just for slightly different file type'
        listnum = [[]]  
        with open(filename) as csvfile:
            lines = csvfile.readlines()
            
            for line in lines[5:]:
                    valuefloat = []
                    values = line.split(',')
                    for value in values:
                        valuefloat.append(float(value))   
                    listnum.append(valuefloat)
                    
                    
            del listnum[0]
                    
        
        
        Voltage = []
        Conductance = []
        Capacitance = []
        InvC = []
        for i in range(len(listnum)):
            
            Voltage.append(listnum[i][0])
            Conductance.append(listnum[i][2])
            Capacitance.append(listnum[i][1])
            InvC.append(1/(listnum[i][2])**2)
            
        
        return Voltage, Capacitance,Conductance,InvC
    
    
    def extractfromallfilenames(self,filenames):
        'Uses extract_parameters function to actually extract across all filenames defined previously'
        Gvalues = [[]]
        Voltagevalues = [[]]
        Capvalues = [[]]
        InvCvalues = [[]]
        for filename in filenames:
            if 'O3' in self.directory: #quick examoke of if statement with other file type
                Voltage,Capacitance, Conductance, InvC = self.extract_parameters_alt(filename)
                Gvalues.append(Conductance)
                Voltagevalues.append(Voltage)
                Capvalues.append(Capacitance)
                InvCvalues.append(InvC)
                
            else:                
                Voltage,Conductance, Capacitance, InvC = self.extract_parameters(filename)
                Gvalues.append(Conductance)
                Voltagevalues.append(Voltage)
                Capvalues.append(Capacitance)
                InvCvalues.append(InvC)
            
            
        del Gvalues[0]
        del Voltagevalues[0]  
        del Capvalues[0]  
        del InvCvalues[0]
        
        
        return Gvalues,Voltagevalues,Capvalues,InvCvalues
        
    

   

    def getpeakGwoverw(self,Gvalues,Capvalues):
        'Inputs experimental parameters in .csv files and uses them to extract values of Dit (density of interface trap density) and Delta_E (energy difference in bandgap from valence/conduction band)'
        

        def Extract(lst):
            'Returns zeroth value of nested list - for getting capacitance and conductance in accumulation only'
            return [item[0] for item in lst]

        C_c = []

        G_c = []
        G_p_w = [[]]
            
        for i in range(len(Gvalues)):  
            G = np.array(Gvalues[i]) #initialise these as values over one frequency at a time to make math more readable
            C = np.array(Capvalues[i])
            freq = np.array(self.omegavalues[i])
            R_s_num =  Extract(Gvalues)[i]**2
            R_s_denom = (Extract(Gvalues)[i])**2 + (freq*(Extract(Capvalues)[i]))**2
            R_s = R_s_num/R_s_denom #calculated series resistance
            C_ox_c = (Extract(Capvalues)[i])*(1 + ((Extract(Gvalues)[i]) / (freq*Extract(Capvalues)[i]))**2) #adjusted oxide capacitance
    
            C_c_num = (G**2 + (freq**2)*(C**2))*C
            C_c_denom = (((G - (G**2 + (freq**2)*(C**2)))*R_s)**2) + (freq**2)*(C**2)
            C_c = C_c_num/C_c_denom  #Adjusted capacitance based on series resistance correction
            
            G_c_num = ((G**2) + (freq**2)*(C**2))*(G - (G**2 + (freq**2)*(C**2))*R_s)
            G_c_denom = ((G - (G**2 + (freq**2)*(C**2))*R_s)**2) + (freq**2)*(C**2)
            G_c = G_c_num/G_c_denom  #Adjusted conductance based on series resistance correction
            
            G_p_w.append((G_c*freq*C_ox_c**2)/(G_c**2 + (freq**2)*(C_ox_c - C_c)**2))
            #Calculate G/omega for each analysed frequency
            
        del G_p_w[0]

        G_p_w_allvoltage = [[]]

        for j in self.voltage_array:
            G_p_w_voltage = []
            for i in range(len(G_p_w)):
            
                G_p_w_voltage.append(G_p_w[i][j]) #Reshape G/omega values to get array conductance values at each voltage point
                
            G_p_w_allvoltage.append(G_p_w_voltage) 
            
            
        del G_p_w_allvoltage[0]
        G_p_w_max = []
        Freq_at_max_G_p_w = []

        for G_p_voltage in G_p_w_allvoltage:
            G_p_w_max.append(max(G_p_voltage)) #Get max conductance value at each voltage point
            Freq_at_max_G_p_w.append(freqvalues_omega[G_p_voltage.index(max(G_p_voltage))]) #Get frequency at which the max G/omega occurs
        
            
        
        v_t = 1e2*(3*k_B*self.temperature/(0.33*m_0))**0.5 #hole thermal velocity
        N = 9.4e14*((self.temperature)**1.5) #effective density of states (done for valence band) cm-3
        tau = 1/np.array(Freq_at_max_G_p_w) #calculate trap lifetime from observed frequency where max G/omega is happening

        Delta_E = k_B*self.temperature*np.log(tau*v_t*sigma*N)/q #Calculate Delta_E (Energy difference in bandgap from valence/conduction band of semiconductor )

        A_cm = A*1e4   #get area in cm^2 instead of m^2
        Dit = np.array(G_p_w_max)*2.5/(q*A_cm) #Calcualte Dit from peak G/omega 


        return Delta_E, Dit
    
    
    
    def returnfromdirectory(self):
        'Returns filepath back to original directory after in preparation for next class to be used with different directory'
       os.chdir('../..')
            
            
    def fitforVfb(self, InvCvalues):
        'linear fit of 1 MHz 1/C^2 vs V curve to extract dopant concentration and flat-band voltage'


        
        InvCap1MHz = InvCvalues[-1] #assumes highest frequency curve is last and is at 1 MHz
        
        
        fit,cov = np.polyfit(Voltagevalues[-1][200:250],InvCap1MHz[200:250],1, cov=True) #linearly fitting about transition from accumulation to depletion on 1 MHz curve may need to change 200:250 if change in where this is located
        
        slope = fit[0]
        intercept = fit[1]
        
        uncertainties = np.sqrt(np.diag(cov))
        
        V_f = np.linspace(-1.7,-0.3,100) #vary these for plotting where transition from depletion to accumulation is
        InvC_f = np.poly1d(fit)
        #plt.plot(V_f,InvC_f(V_f),color = 'r',linestyle = '--') #can plot if plot of 1/C^2 vs V is there as well
        
        
        
        
        V_fb = -intercept/slope - k_B*self.temperature/q #Flat-band voltage occurs at linear extrapolation to 1/C^2 = 0
        
        N_d = 2/(epsilon_Ge*epsilon_0*q*(A**2)*slope*1e6) #Dopant concentration in cm-3
        
        return V_fb, N_d



#%%


#Example class definition

values77K = Ditextract('S40654_GeCAP_GeO2_SiN_FG380C1hr/77K','*_GeCAP_*_*Hz_500um.csv', 77,freqvalues_omega,Voltage_index_77)



#%%

def getparameters(obj):
    'Extracts raw experimental data - Voltage, Conductance, Capacitance and InvC (1/C^2)'
    obj.gotodirectory()
    filenames = obj.extractfilenames(obj.filenames_extracted)
    Gvalues,Voltagevalues,Capvalues,InvCvalues = obj.extractfromallfilenames(filenames)
    
    obj.returnfromdirectory()

    return Voltagevalues, Gvalues, Capvalues, InvCvalues



def getDitDeltaE(obj):
    'Defines order in which class functions are carried out to get Dit and Delta_E'
    obj.gotodirectory()
    filenames = obj.extractfilenames(obj.filenames_extracted)
    Gvalues,Voltagevalues,Capvalues,InvCvalues = obj.extractfromallfilenames(filenames)
    Delta_E, Dit = obj.getpeakGwoverw(Gvalues, Capvalues)
    obj.returnfromdirectory()
    
    return Delta_E, Dit



    
    
#Incidence of functions
Voltagevalues, Gvalues, Capvalues, InvCvalues = getparameters(values77K)
Delta_E, Dit = getDitDeltaE(values77K)


#plot figure of Dit vs Delta_E
plt.figure()
plt.scatter(Delta_E, Dit,color='tab:blue')
plt.xlabel('Energy (eV)')
plt.ylabel('D$_{it}$ (eV$^{-1}$ $cm^{-2}$) ')
plt.title('D$_{it}$ vs Energy')
ax = plt.gca()
ax.set_yscale('log')
plt.grid(which = 'major',linestyle = '-')
plt.grid(which = 'minor',linestyle = 'dotted')





