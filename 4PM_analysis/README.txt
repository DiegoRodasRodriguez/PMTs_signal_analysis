****************************************************************
****************************************************************
**********************                    **********************
**********************    4PM_analysis    **********************
**********************                    **********************
****************************************************************
****************************************************************


	-----> Specs to run the main scripts: <-----


·MATLAB version			: 	v9.14.0.2254940 (R2023a) Update 2
					v9.8.


·Add-ons needed			: 	Optimization Toolbox
					SimBiology


					* They could be not necessary (maybe not both of them) for this set of scripts but they are 
					mandatory for other scripts, so I'd recommend to install them.


·Internal paths			:	In the description of every script there will be a section 'Paths to change' specifying 
					the line and var name which are need to change.
		

(__***__)								(__***__)
(__***__)	A LOT OF WARNINGS APPEAR WHEN RUNNING THIS SCRIPTS	(__***__)
(__***__)								(__***__)


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------ MAIN SCRIPTS ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------

2 initialization scripts (AnaSwanArCF4_S1ini, AnaSwanArCF4_S2ini) + 2 (main) analysis scripts (AnaSwanArCF4_S1_PRO, AnaSwanArCF4_S2_PRO) . 


There are 2 scripts from each kind due to the 2 types of signals we want to analyze S1 and S2.

There are data-files with S1 only and data-files that contains both S1 and S2 so they allow to estimate the drift velocity.



·AnaSwanArCF4_S1ini:

	
	Sets the values of some parameters that include info about the experimental measures or how they had been taken.		
				(suppression of low S/N levels, QE for Ar or Xe, corrections,...)

		-With the lower StoN suppression we are performing a cut, so it could occur that this default values are not general 
		for every data-set.


	It also specifies where the data files are saved with the DIR path.			(DIR = [loc_path, general_path])

		-It is necessary to change this path to your local disk path of Dropbox (loc_path) in every single path 
		('C:\Users\user_name\RareEventsGroup Dropbox\...'). To make this more comfortable DIR was divided into a local
		 path that would depend on the computer where the script will be running and a general path that is the same 
		for everyone in the Dropbox group:
		
			-loc_path: 'C:\Users\user_name\RareEventsGroup Dropbox\'
			-general_path: '...\HOME_RareEventsGroup\...\...'


	This script also contains info about the CALIBRATION of the PMs in terms of accumulated charge to photons equiv. (vars: QPM1toPh...)




	---Paths to change---
		
		line 40		:	loc_path
 	
	   
·AnaSwanArCF4_S2ini:
	
	Very similat to the previous one. It sets some parameters and makes StoN suppression and also specify the loc_path again 
	for the data-sets.
	
		- There are 2 parameters that there are not in the S1 one to perform a relative normalization to estimate 
		position and other to shift all avg wvfs to match PM1. Both of this are '= 1' by default.
	

	---Paths to change---
		
		line 44		:	loc_path


·AnaSwanArCF4_S1_PRO:

	It performs all the analysis. The Data and Workflow are very good described in the NOTES section at the beggining of 
	the script.

	
	Input	:	Data 
	
	Output	:	Summary + 21 plots.


					(MAKE HERE A SUMMARY ABOUT THE PERFORMED ANALYSIS)**

·AnaSwanArCF4_S2_PRO:

	Almost the same that S1. (Briefly comment about the differences when I take a deep look into the script)
	

	Input	:	Data 
	
	Output	:	Summary + 19 plots.

	

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




**********************************************************************************************************************************************
**********************************************************************************************************************************************
**********************************************************************************************************************************************
**********************************************************************************************************************************************
**********************************************************************************************************************************************
**********************************************************************************************************************************************
THIS SECTION IS TEMPORARY....... :)

--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------  TO DO NOTES  -----------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------


There a lot of TO DO in the analysis scripts. I will be adding them into here when I 
start to check it more carefully. 


	-----> TO DELETE <-----


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------   DOUBTS     ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------


--> Why is the trigger referring to PM1 and PM4 ???????
--> 

