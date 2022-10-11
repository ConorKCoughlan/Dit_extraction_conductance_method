# Dit_extraction_conductance_method

This script can be used to calculate the density of interface trap densities (Dit) as a function of energy within the bandgap of a semiconductor material (i.e. the value ΔE where ΔE =  E - Ev or Ec - E) using the conductance method. For background reading of this method, please see the following references:

 https://doi.org/10.1063/1.3520431
 
 MOS (Metal Oxide Semiconductors): Physics and Technology -  E H Nicollian; J R Brews - Free Download Here (https://vdoc.pub/documents/mos-metal-oxide-semiconductor-physics-and-technology-669oh75ou100)
 
 
## Note on raw data format
 This script starts from extracting raw data from Capacitance-Voltage (CV) and Conductance Voltage (GV) measurements while both the AC electrical frequency and the temperature are varied. The raw data was measured on a cryogenically cooled Keysight B1500a Semiconductor Parameter Analyser, and the raw data format reflects this. On a different system, it is likely this would need to be adjusted to account for different data format. 
 
 ## How to use the code
 The class Dit extract contains all the functions required to extract Dit from raw data, including navigating to the relevant directories, extracting and organising the raw data and then calculating the Dit and Delta_E parameters subsequently. Class inputs can be altered by changing any of the input parameters into the class (Directory, filenames_extracted, temperature, frequency, analysed voltage points). 
 
 e.g. values77K = Ditextract('S40654_GeCAP_GeO2_SiN_FG380C1hr/77K','*_GeCAP_*_*Hz_500um.csv', 77,freqvalues_omega,Voltage_index_77)
 
 ### Frequency
 As often a range of frequency values are taken for every temperature point, the parameter 'freqvalues_omega' accounts for all inputted frequency at the respective temperature and may be modified directly (line 13). 
 
 
 ### Voltage
 Similarly for the analysed voltage points parameters (given as Voltage_index_temperature (lines 19-21), they can be modified here rather than inputting directly into class. This is crucial as only certain voltage points will contain a maximum G/omega vs omega conductance peak within the middle of frequency range, and it is only those analysed voltage values that will determine the Dit, as the conductance peak is occurring inside the measured frequency range. This often requires a little bit of trial and error with selecting different ranges of voltage. 
 
 
 ## Additional functions
 
 This code also includes the functions getparameters(obj) and getDitDeltaE(obj), which will take the class and extract the raw parameters and Dit and DeltaE respectively. Also included within the class, there is the function fitforVfb, which will take the high frequency (1 MHz) 1/C^2 vs V plot and extract the V_fb (flat-band voltage) and N_d (dopant concentration) of the curve using a linear fit on the transition from accumulation to depletion. This analysis is typically done at room temperature. 
 
 
 
 
 
 


