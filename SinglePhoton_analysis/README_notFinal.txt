**********************************************************************
**********************************************************************
**********************                          **********************
**********************  SINGLE PHOTON ANALYSIS  **********************
**********************                          **********************
**********************************************************************
**********************************************************************



Specs to run the main scripts:

·MATLAB version: 	v9.14.0.2254940 (R2023a) Update 2
			v.9.8.

·Add-ons needed: 	Optimization Toolbox


·Internal paths			:	In the description of every script there will be a section 'Paths to change' specifying 
					the line and var name which are need to change.



------------------------------------------
------------- MAIN FUNCTIONS -------------
------------------------------------------

2 scripts + 1 main script: SinglePhotonAna_PRO, SinglePhotonIni,ReadDataCAEN_PRO



·SinglePhotonIni:
	
	-Define DIR and FILE where data to be read is stored.
	-Add paths where the aux-fnc are stored.
	-Define the function used to read data which depends of the hardware (DataType): 'CAEN' or 'OSCI'.
		'CAEN' --- Data taken with CAEN card
		'OSCI' --- Data taken with Tektronix scope
	-Define some parameters and windows to look for Amplitude and Charge signals and nChannels (nCh).


	---Paths to change---
		
		line 302		:	loc_path


·ReadDataCAEN_PRO:

	Function used to read data which has been taken with the CAEN card.


·Main:	SinglePhotonAna_PRO 

	It does all the analysis things. I will be more speciffic about this in other .txt.



--------------------------------------------------------------

Auxiliar scripts that have been used are called at the beggining each analysis script which the 'addpath' func 
referring their corresponding path in your local directory or Dropbox one.




------------------------------------------
------------------- DATA -----------------
------------------------------------------

The data is saved in the 'Data' folder. 
This dataset was recorded with the CAEN card.








