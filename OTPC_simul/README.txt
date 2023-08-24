***************************************************************
***************************************************************
**********************                   **********************
**********************  OTPC SIMULATION  **********************
**********************                   **********************
***************************************************************
***************************************************************



	-----> Specs to run the main scripts: <-----

·MATLAB version			: 	v9.14.0.2254940 (R2023a) Update 2
					v9.8.


·Add-ons needed			: 	SimBiology*
					Statistics and Machine Learning Toolbox (I think it is included when you install SimBiology)


·Internal paths			:	In the description of every script there will be a section 'Paths to change' specifying 
					the line and var name which are need to change.


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------ MAIN SCRIPTS ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------

6 scripts + 1 main script (simulation): 	CMOSresponse, GenStraightTrack, OTPC_simul, PMTresponse, TPCanode, TPCdrift, TPCprimaryResponse

It is not necessary to run all the scrips before the main one. When you run the main script it will call the others so if any of them fails it will print a console-log with the error (Script + line).	
	
	
		
					_____DEFINITION_____											_____VARIABLES_____



·TPCprimaryResponse:

	Simulates primary response from TPC generated for the first ionization and later scintillation. 

	Input (5)		:	position and time of the energy deposites.								(x, y, z, t, DeltaE)	

	Output (9)		: 	position and time of the e- generated from 1st ionization and time and wavelenght from scintillation.	(xe, ye, ze, te, xph, yph, zph, tph, Lph)

	Params/global (7)	:	specs needed for ionization and scintillation (Fano factors and more) 					(Wi, Ws, tau1, tau2, Ratio12, FanoQ, FanoS)


·TPCdrift:

	Simulates the drift suffered by the first e- in the TPC.

	Input (4)		: 	position and time of the initial e- (generated from 1sr ionization // TPCprimaryResponse output).	(xe, ye, ze, te)

	Output (4)		: 	position and time of the e- at Anode (input for TPCanode).		 				(xeA, yeA, zeA, teA)

	Params/global (4)	:	drift/diffussion specs. It depends of the gas used (?)							(P, vd, DL, DT)


·TPCanode:

	Simulates the Anode response. Anode is where the amplification system is located (GEMs) and electrons are being multiplicated.

	Input (4)		: 	position and time of the photons produced at Anode. 							(xeA, yeA, zeA, teA)

	Output (5)		: 	position and time of the e- produced in the Anode and statistics weights of those ph. 			(xePh, yePh, zePh, tePh, WePh)

	Params/global (8)	:	specs from the amplification system used (GEM, mesh...)							(hole, pitch, gap, OptGain, vdGap, tau1Gap. tau2Gap, Ratio12Gap)


·GenStraightTrack:

	Simulates straight tracks in a given material. Material and particle whose track will be simulated info are given in the inputs as M, Range, dEdx,... 
	
	Inputs (7)		:	Initial position (cilindric) and interaction particle-material 						(z0, theta0, phi0, Range, dEdx, M, nsteps)

	Output (5)		:	Final track (position and time) and energy deposit. 							(x, y, z, t, DeltaE)


·CMOSresponse:
	
	Simulates the response of a CMOS sensor as it records the accumulated charge.
	
	Inputs (5)		:	Position, time and weight for S2 									(xePh, yePh, zePh, tePh, WePh)

	Output (1)		:	X-Y matrix with nphotons as entries. 									(NphXY) 

	Params/global (7+1)	:	Speciffications of the camera and optical details 							(M, lensN, sigmaNph, QECMOS, rebin, Npixel, pixelSize, +Tmesh)



·PMTresponse:

	Simulates the PMTs response to a photons interaction. It will simulate both, S1 and S2 signals expected.

	Input (9)		:	Info about the S1 and S2 signals (position, time and weight). 						(S2: xS2, yS2, zS2, tS2, WS2  ;  S1: xS1, yS1, zS1, tS1)

	Output (12)		:	PMs waveforms (wvf) and times for every single PM (4: up, down, left, right) for S1 and S2 signals. 	(PML_wvf, PMR_wvf, PMU_wvf, PMD_wvf, Twvf, tPML_S1, tPMR_S1,tPMU_S1, tPMD_S1, tPML_S2, 																				tPMR_S2, tPMU_S2, tPMD_S2)

	Params/global (12+2)	:	position and specs of the PMs. 										(zPM, rPM, RPM, QEPM, Tmesh, sigmaAovA1, sigmaNoise, A1, sigmaTresponse, meanTresponse, Tbin, 																		BW, +gap, +vdGap)



-----------------------------------------

·Main:			OTPC_simul

	It calls previous functions and defines the params and specs needed. It also provide the code lines to plot the info and show us 6 plots and a Resume.
	
	Input			:	Previous functions whose has been described before.
	
	Output			:	Resume and 6 plots 											(primary ionization (PI), primary scintillation (PS), PI@anode, PS@anode, PM signals, CMOS 																			capture)

	---Paths to change---
		
		line 39		:	loc_path

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

No data is needed to run this script.



***************************************************************************************************************************************************************************************************************************************
***************************************************************************************************************************************************************************************************************************************
***************************************************************************************************************************************************************************************************************************************
***************************************************************************************************************************************************************************************************************************************
***************************************************************************************************************************************************************************************************************************************
***************************************************************************************************************************************************************************************************************************************

--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------  TO DO NOTES  -----------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------

This is only a copy&paste from the scripts notes:



· TPCprimaryResponse:
	
	-Include scintillation spectrum
	-Fano included for constant energy loss (modify when B-B included!)	

· TPCdrift:

	-Include initialization parameters from file
	-Add attachment

· TPCanode:
	
	-It assumes Furry law to calculate statistical weight, assigned to anode.
	-Scintillation during e- transit included through a convolution in PMTresponse (large numerical overhead otherwise)

· GenStraightTrack:

	-Include Bethe-Bloch

· CMOSresponse:
	
	-Include/introduce function for weighted histogram to include WePh electron by electron.
	-In this code also exists a temporary patch!!! Look for improve this.

· PMTresponse:

	-Gap transit time introduced through a convolution with a flat distribution. True for EL. For CF4 an exponential will be more suited. Include this as a parameter.
	-No reflection included, only direct light.

	REFERENCE SYSTEM:

									      ^ y
		       o	(UP)				     	      |
	      (LEFT) o   o 	(RIGHT)					x <---|
		       o	(DOWN)

· OTPC_simul:

	-Include a more realistic energy-loss
	-PMT time response totally approximate
	-SinglePhoton width and noise from Q-analysis, not from A-analysis.
	-Negative Tail of SinglePhoton probably unrealistic.
	-Include PMT-by-PMT variations.
	-Include global response parameters as workspaces and document, in order to understand the various assumptions. All parameters should be revisited then.
	-Gap transit-time included as a convolution with a flat response. Allow to change between exponential or flat, depending on the case (CF4, noble).
	-No CMOS offset. It requires maps...
	-CMOS noise approx through normal distribution.

	15s with rebin 10, 10 min with rebin 1. -----> Time elapsed.

	Re-check CMOS numbers.



	-----> TO DELETE <-----


--------------------------------------------------------------
--------------------------------------------------------------
------------------------ 	      ------------------------
------------------------   DOUBTS     ------------------------
------------------------ 	      ------------------------
--------------------------------------------------------------
--------------------------------------------------------------


--> Why only 1 ph per e-? Is it a suppose??
--> Furry law?? -- stat weigh for anode
--> S2 weight ?? 
--> Output CMOS
--> TPCdrift: drif and diffusion depends on the gas used??


