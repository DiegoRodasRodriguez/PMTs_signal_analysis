****************************************************************
****************************************************************
**********************                    **********************
**********************    1PM_analysis    **********************
**********************                    **********************
****************************************************************
****************************************************************



*** Pay special attention to the following symbol (__***__)



	-----> Specs to run the main scripts: <-----


·MATLAB version			: 	v9.14.0.2254940 (R2023a) Update 2
					v9.8.
					
					***some versions may cause problems mainly with the ReadDataCAEN_PRO script while 
					loading data


·Add-ons needed			: 	Optimization Toolbox
					SimBiology*

					* They could be not necessary (maybe not both of them) for this set of scripts but 
					they are mandatory for other scripts, so I'd recommend to install them.



·Internal paths			:	In the description of every script there will be a section 'Paths to change' specifying 
					the line and var name which are need to change.



--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------ MAIN SCRIPTS ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------


3 main (FATGEMsetrun_PRO, FATGEMini_PRO, FATGEMana_PRO) + 3 data readout scripts (ReadDataCAEN*, ReadDataCAEN_PRO, ReadDataOsciSingle) ----->(6 total)

(__***__) MAIN SCRIPT: 	FATGEMana_PRO		(this one calls the others, so this is the only one you need to run.)


*I think this one is not mandatory. It's some kind of an old version of the PRO so it should be deleted in future versions.



·FATGEMsetrun_PRO:

	Defines the directory 'DIR' where the data to use in the later analysis is saved.
	
	
	---Paths to change---
		
		line 125		:	loc_path



·FATGEMini_PRO:

	Defines some statistics values needed for the analysis (cuts, comparisons, ...) whose depends of the loaded dataset.

	Basically it compares the DIR with all of the datasets saved in Dropbox and define some values depending on the DIR.
	


·FATGEMana_PRO:

	Performs all the analysis. It calls the previous scripts in order to get the data, set the statistical params and 
	permforms the analysis.

	(__***__) 	To understand the analysis (check _____ file). -----> (It will be explained in another .txt file.)

	

	Input 	:	Data

	Output	:	Summary with relevant data (energy resolutions, n events, charge, witdths,...) and 11 plots.
				
			Plots: 
			1-> PMT wvfs, 		
			2-> pulse info, 	
			3 to 7-> charge info, 		
			8 -> drift,
			9-> (?), 	
			10-> energy info, 	
			11-> (?)



------------------------ Data readout ------------------------ 


·ReadDataCAEN_PRO:

	Script to read data taken with the CAEN card.


·ReadDataOSciSingle:

	Script to read data taken with the Osciloscope.


·ReadDataCAEN *(It will be deleted): 

	Script to read data taken with the CAEN card. I think this one is not necessary (maybe an old version of the PRO).


--------------------------------------------------------------

Auxiliar scripts that have been used are called at the beggining each analysis script which the 'addpath' func 
referring their corresponding path in your local directory or Dropbox one.



--------------------------------------------------------------
--------------------------------------------------------------
---------------------------- 	  ----------------------------
---------------------------- DATA ----------------------------
---------------------------- 	  ----------------------------
--------------------------------------------------------------
--------------------------------------------------------------


Data used to test this scripts are stored in the group Dropbox.

If you want to change the data you only to add the corresponding path in the initialization scripts.



*****************************************************************************************************************************
*****************************************************************************************************************************
*****************************************************************************************************************************
*****************************************************************************************************************************
*****************************************************************************************************************************
*****************************************************************************************************************************


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------  TO DO NOTES  -----------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------

· Test the current scripts (with the 'eval' removal) in the MATLAB v9.8.   (the current ones has been only tested in the v9.14)
· Complete the PLOTS info


	-----> TO DELETE <-----


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------   DOUBTS     ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------


--> 
--> 
--> 
--> 
--> 

