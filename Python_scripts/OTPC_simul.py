# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 11:51:38 2023

@author: diego
"""


"""
We are going to rewrite the MATLAB code into Python and using OOP with classes.

There are 7 diff matlab scripts that there will be classes in this new 
redit code.

Main class is going to be OTPC_simul which is the one 

"""

# Import the modules

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

pi=math.pi


def hist1D(y,x):
    
    N, _ =np.histogram(y, bins=x)
    
    if np.size(x) != np.size(N):
        N=N.T
    
    plt.stairs(N,x)
    plt.show
        
    return N

"""
def PyConvolve(U1,X1,U2,X2, Mode, Type='FFT', verbose=1):
    
    sizeU1=np.size(U1)
    sizeX1=np.size(X1)
    sizeU2=np.size(U2)
    sizeX2=np.size(X2)
    
    if sizeU1[0] != 1:
        U1=U1.T
        
    if sizeX1[0] != 1:
        X1=X1.T
        
    if sizeU2[0] != 1:
        U2=U2.T
        
    if sizeX2[0] != 1:
        X2=X2.T
    
    return I, Time
"""

#def hist3D():
    

class OTPCsimulator:
    
    
    def GenStraightTrack(z0, theta0,phi0,Range,dEdx, M,nsteps):
        
        #This returns x,y,z,t,DeltaE in matlab
        
        
        Ekin=dEdx*Range         # [MeV]
        clight=30               # [cm/ns]
        
        rlength=Range*np.random.rand(nsteps)
        rlength.sort()
        beta=np.sqrt(2*Ekin/M)
        vel=beta*clight
        
        t=np.zeros(len(rlength))
        totLength=0
        
        for i in range(2,len(rlength)):
            
            deltaLength=rlength[i]-rlength[i-1]
            totLength=totLength+deltaLength
            t[i]=t[i-1]+deltaLength/vel
            beta=np.sqrt(2*(Ekin-dEdx*totLength)/M)
            vel=beta*clight
            
        x=rlength*np.sin(theta0*pi/180)*np.cos(phi0*pi/180)
        y=rlength*np.sin(theta0*pi/180)*np.sin(phi0*pi/180)
        z=z0+rlength*np.cos(theta0*pi/180)
        DeltaE=Ekin/nsteps*np.ones(len(rlength))*(1e6)
        
        return x,y,z,t,DeltaE                                                   # THIS OUTPUT COULD BE A LIST [POS,TIME] AND DELTA E
    

    def TPCprimaryResponse(x,y,z, t, DeltaE):
        
        global Wi, Ws, tau1, tau2, Ratio12, FanoQ, FanoS
        
        neToTmean=round(sum(DeltaE)/Wi)
        nphToTmean=round(sum(DeltaE)/Ws)
        
        sigmaStepQ=(np.sqrt(FanoQ)*np.sqrt(neToTmean))/np.sqrt(len(DeltaE))
        sigmaStepS=(np.sqrt(FanoS)*np.sqrt(nphToTmean))/np.sqrt(len(DeltaE))
        
        ne=np.zeros(len(DeltaE))
        nph=np.zeros(len(DeltaE))
        
        neToT=0
        nphToT=0
        
        for i in range(len(DeltaE)):
            
            ne_mean=DeltaE[i]/Wi
            ne[i]=round(np.random.normal(ne_mean,sigmaStepQ))
            nph_mean=DeltaE[i]/Ws
            nph[i]=round(np.random.normal(nph_mean,sigmaStepS))
            neToT=neToT+ne[i]
            nphToT=nphToT+nph[i]
        
        xe=np.zeros(int(neToT))
        ye=np.zeros(int(neToT))
        ze=np.zeros(int(neToT))
        te=np.zeros(int(neToT))
        
        xph=np.zeros(int(nphToT))
        yph=np.zeros(int(nphToT))
        zph=np.zeros(int(nphToT))
        tph=np.zeros(int(nphToT))
        
        Lph=np.zeros(int(nphToT))
        
        I0=0
        nphCumm=0
        
        
        #Loop over energy deposits
        for i in range(len(DeltaE)):
            
            #electrons
            If=I0+ne[i]
            xe[int(I0):int(If)]=x[i]
            ye[int(I0):int(If)]=y[i]
            ze[int(I0):int(If)]=z[i]
            te[int(I0):int(If)]=t[i]
            I0=If+1
            
            #photons
            for j in range(int(nph[i])):
                
                rand_aux=np.random.rand()
                
                if rand_aux > Ratio12:
                    
                    tph[int(nphCumm)+j]=np.random.exponential(tau2)
                    
                else:
                    tph[int(nphCumm)+j]=np.random.exponential(tau1)
                
                
                xph[int(nphCumm)+j]=x[i]
                yph[int(nphCumm)+j]=y[i]
                zph[int(nphCumm)+j]=z[i]
                Lph[int(nphCumm)+j]=-1
                
            nphCumm=nphCumm+nph[i]
            
        return xe, ye, ze, te, xph, yph, zph, tph, Lph
    
    
    def TPCdrift(xe,ye,ze,te):
        
        global P, vd, DL, DT
        
        xeA=np.zeros(np.size(xe))
        yeA=np.zeros(np.size(ye))
        zeA=np.zeros(np.size(ze))
        teA=np.zeros(np.size(te))
        
        for i in range(len(xe)):
            
            teA[i]=te[i]+np.random.normal(loc=ze[i]/vd,scale=(DL*np.sqrt(ze[i]))/vd)
            
            xeA[i]=np.random.normal(loc=xe[i],scale=DT*np.sqrt(ze[i]))
            
            yeA[i]=np.random.normal(loc=ye[i],scale=DT*np.sqrt(ze[i]))
            
            zeA[i]=0
        
        return xeA,yeA,zeA,teA
                    
    
    def TPCanode(xeA,yeA,zeA,teA):
        
        global hole, pitch, gap, OptGain, vdGap, tau1Gap, tau2Gap, Ratio12Gap
        
        xePh=pitch*np.rint(xeA/pitch)
        yePh=pitch*np.rint(yeA/pitch)
        zePh=zeA
        
        tePh=np.zeros(np.shape(zeA))
        WePh=zePh
        
        for i in range(len(xeA)):
            
            #smear within the hole
        
            r2_rnd=np.random.rand()
            r_rnd=np.sqrt(r2_rnd)*hole/2
            phi_rnd=2*pi*np.random.rand()
            xePh[i]=xePh[i]+r_rnd*np.sin(phi_rnd)
            yePh[i]=yePh[i]+r_rnd*np.cos(phi_rnd)
            #Add avalanche statistics
            WePh[i]=np.random.exponential(OptGain)
            
            if np.random.rand() > Ratio12Gap:
                
                tePh[i]=teA[i]+np.random.exponential(tau2Gap)
            
            else:
                
                tePh[i]=teA[i]+np.random.exponential(tau1Gap)
        
        
        return xePh, yePh, zePh, tePh, WePh
    
    
    def CMOSresponse(xePh,yePh,WePh):
        
        global M, lensN, sigmaNph, QECMOS, rebin, Npixel, pixelSize, Tmesh
        
        #Calculate Solid angle and sensor size
        
        OmegaCMOS=1/(16*lensN**2)*1/(1/M+1)**2
        WCMOS=Tmesh*OmegaCMOS*QECMOS
        
        xsideCMOS=Npixel*pixelSize
        ysideCMOS=Npixel*pixelSize
        
        #Project anode positions
        xeCMOS=M*xePh   #Inverts!
        yeCMOS=M*yePh   #Inverts!
        
        #Create matrix of weights
        Xmax=xsideCMOS/2-rebin*pixelSize
        Ymax=ysideCMOS/2-rebin*pixelSize    #1 bin needs to be substracted in order to give the right cound
        
        Xmin=-xsideCMOS/2
        Ymin=-ysideCMOS/2
        
        x_edg=np.arange(Xmin,Xmax,rebin*pixelSize)
        y_edg=np.arange(Ymin,Ymax,rebin*pixelSize)
        NphXY,xedge,yedges=np.histogram2d(xeCMOS, yeCMOS, bins=(x_edg,y_edg))
        
        sizeNphXY=np.shape(NphXY)
        NpixXrebin=sizeNphXY[0]
        NpixYrebin=sizeNphXY[1]
        
        #Include statistical weights and noise
        
        WePhMean=np.mean(WePh)     #Note: temporary patch
        
        for i in range(NpixXrebin):
            for j in range(NpixYrebin):
                
                NphXY[i,j]=np.random.poisson(WePhMean*WCMOS*NphXY[i,j])+np.random.normal(0,sigmaNph*np.sqrt(rebin*rebin))
        
        
        return NphXY
    
    
    
    def PMTresponse(xS2, yS2,zS2,tS2, WS2, xS1, yS1, zS1, tS1):
        
        global zPM, rPM, RPM, QEPM, Tmesh, sigmaAovA1, sigmaNoise, A1, sigmaTresponse, meanTresponse, Tbin, BW, gap, vdGap
        
        clight=30   #[cm/ns]
        
        #######################################################################
        ############   GENERATE TIME ENTRIES FOR DETECTED PHOTONS   ###########
        #######################################################################
        
        #Primary scintillation weights/probabilities per time entry
        WS1L=np.zeros(np.shape(tS1))
        WS1R=np.zeros(np.shape(tS1))
        WS1U=np.zeros(np.shape(tS1))
        WS1D=np.zeros(np.shape(tS1))
        
        lengthS1L=np.zeros(np.shape(tS1))
        lengthS1R=np.zeros(np.shape(tS1))
        lengthS1U=np.zeros(np.shape(tS1))
        lengthS1D=np.zeros(np.shape(tS1))
        
        for i in range(len(tS1)):
            
            lengthS1L[i]=np.sqrt((zS1[i]-zPM)**2+(xS1[i]-RPM)**2+yS1[i]**2)
            lengthS1R[i]=np.sqrt((zS1[i]-zPM)**2+(xS1[i]-RPM)**2+yS1[i]**2)
            lengthS1U[i]=np.sqrt((zS1[i]-zPM)**2+xS1[i]**2+(yS1[i]-RPM)**2)
            lengthS1D[i]=np.sqrt((zS1[i]-zPM)**2+xS1[i]**2+(yS1[i]+RPM)**2)
            
            WS1L[i]=(1/4)*rPM**2/lengthS1L[i]**2
            WS1R[i]=(1/4)*rPM**2/lengthS1R[i]**2
            WS1U[i]=(1/4)*rPM**2/lengthS1U[i]**2
            WS1D[i]=(1/4)*rPM**2/lengthS1D[i]**2
            
        WS1L=WS1L*QEPM*Tmesh
        WS1R=WS1R*QEPM*Tmesh
        WS1U=WS1U*QEPM*Tmesh
        WS1D=WS1D*QEPM*Tmesh
            
        
        #Secondary scintillation weights/probabilities per time entry
        WS2L=np.zeros(np.shape(tS2))
        WS2R=np.zeros(np.shape(tS2))
        WS2U=np.zeros(np.shape(tS2))
        WS2D=np.zeros(np.shape(tS2))
        
        lengthS2L=np.zeros(np.shape(tS2))
        lengthS2R=np.zeros(np.shape(tS2))
        lengthS2U=np.zeros(np.shape(tS2))
        lengthS2D=np.zeros(np.shape(tS2))
        
        for i in range(len(tS2)):
            
            lengthS2L[i]=np.sqrt((zS2[i]-zPM)**2+(xS2[i]-RPM)**2+yS2[i]**2)
            lengthS2R[i]=np.sqrt((zS2[i]-zPM)**2+(xS2[i]-RPM)**2+yS2[i]**2)
            lengthS2U[i]=np.sqrt((zS2[i]-zPM)**2+xS2[i]**2+(yS2[i]-RPM)**2)
            lengthS2D[i]=np.sqrt((zS2[i]-zPM)**2+xS2[i]**2+(yS2[i]+RPM)**2)
            
            WS2L[i]=(1/4)*rPM**2/lengthS2L[i]**2
            WS2R[i]=(1/4)*rPM**2/lengthS2R[i]**2
            WS2U[i]=(1/4)*rPM**2/lengthS2U[i]**2
            WS2D[i]=(1/4)*rPM**2/lengthS2D[i]**2
            
        WS2L=WS2L*QEPM*Tmesh
        WS2R=WS2R*QEPM*Tmesh
        WS2U=WS2U*QEPM*Tmesh
        WS2D=WS2D*QEPM*Tmesh
        
        
        """
        NOTE: a x-check in the number of photons is here in order to avoid 
        mistakes in the statistics. Given the very low weight (much less than 1 
        since each time represents an individual photon or electron), the 
        binomial approximation will be generally fine.
        Note that PML will detect slightly more photons through this algorithm 
        (it is interrogated first), and will only detect one photon per weighted 
        photon. These approximations are not expected to have practical consequ
        ences in general.
        
        """
        
        kL_S1=0
        kR_S1=0
        kU_S1=0
        kD_S1=0
        
        kL_S2=0
        kR_S2=0
        kU_S2=0
        kD_S2=0
        
        tPML_S1=np.zeros(len(tS1))
        tPMR_S1=np.zeros(len(tS1))
        tPMU_S1=np.zeros(len(tS1))
        tPMD_S1=np.zeros(len(tS1))
        
        tPML_S2=np.zeros(len(tS2))
        tPMR_S2=np.zeros(len(tS2))
        tPMU_S2=np.zeros(len(tS2))
        tPMD_S2=np.zeros(len(tS2))
        
        tPML_S1[0]=-99999
        tPMR_S1[0]=-99999
        tPMU_S1[0]=-99999
        tPMD_S1[0]=-99999
        
        tPML_S2[0]=-99999
        tPMR_S2[0]=-99999
        tPMU_S2[0]=-99999
        tPMD_S2[0]=-99999
        
        
        for i in range(len(tS1)):
            
            if np.random.rand() <= WS1L[i]:
                tPML_S1[kL_S1+1]=tS1[i]+lengthS1L[i]/clight
                kL_S1=kL_S1+1
            
            elif np.random.rand() <= WS1R[i]:
                tPMR_S1[kR_S1+1]=tS1[i]+lengthS1R[i]/clight
                kR_S1=kR_S1+1
            
            elif np.random.rand() <= WS1U[i]:
                tPMU_S1[kU_S1+1]=tS1[i]+lengthS1U[i]/clight
                kU_S1=kU_S1+1
            
            elif np.random.rand() <= WS1D[i]:
                tPMD_S1[kD_S1+1]=tS1[i]+lengthS1D[i]/clight
                kD_S1=kD_S1+1
        
        
        for i in range(len(tS2)):
            
            if np.random.rand() <= WS2L[i]:
                tPML_S2[kL_S2+1]=tS2[i]+lengthS2L[i]/clight
                kL_S2=kL_S2+1
            
            elif np.random.rand() <= WS2R[i]:
                tPMR_S2[kR_S2+1]=tS2[i]+lengthS2R[i]/clight
                kR_S2=kR_S2+1
            
            elif np.random.rand() <= WS2U[i]:
                tPMU_S2[kU_S2+1]=tS2[i]+lengthS2U[i]/clight
                kU_S2=kU_S2+1
            
            elif np.random.rand() <= WS2D[i]:
                tPMD_S2[kD_S2+1]=tS2[i]+lengthS2D[i]/clight
                kD_S2=kD_S2+1
                
                
        #######################################################################
        ########################   GENERATE WAVEFORMS   #######################
        #######################################################################
        
        #Generate waveforms (I. Add time entries)
        Twvf=np.arange(-1000,1000,Tbin)
        
        PML_wvf_S1=np.zeros(np.shape(Twvf))
        PMR_wvf_S1=np.zeros(np.shape(Twvf))
        PMU_wvf_S1=np.zeros(np.shape(Twvf))
        PMD_wvf_S1=np.zeros(np.shape(Twvf))
        
        PML_wvf_S2=np.zeros(np.shape(Twvf))
        PMR_wvf_S2=np.zeros(np.shape(Twvf))
        PMU_wvf_S2=np.zeros(np.shape(Twvf))
        PMD_wvf_S2=np.zeros(np.shape(Twvf))
        
        
        if kL_S1 > 0:
            PML_wvf_S1=hist1D(tPML_S1, Twvf)
        if kR_S1 > 0:
            PMR_wvf_S1=hist1D(tPMR_S1, Twvf)
        if kU_S1 > 0:
            PMU_wvf_S1=hist1D(tPMU_S1, Twvf)
        if kD_S1 > 0:
            PMD_wvf_S1=hist1D(tPMD_S1, Twvf)
            
            
        if kL_S2 > 0:
            PML_wvf_S2=hist1D(tPML_S2, Twvf)
        if kR_S2 > 0:
            PMR_wvf_S2=hist1D(tPMR_S2, Twvf)
        if kU_S2 > 0:
            PMU_wvf_S2=hist1D(tPMU_S2, Twvf)
        if kD_S2 > 0:
            PMD_wvf_S2=hist1D(tPMD_S2, Twvf)
        
        
        #Generate wavefrorms (II. Add PM statistics)
        
        for i in range(len(Twvf)-1):
            
            if PML_wvf_S1[i] > 0:
                PML_wvf_S1[i]=np.random.normal(PML_wvf_S1[i], sigmaAovA1*np.sqrt(PML_wvf_S1[i]))
            if PMR_wvf_S1[i] > 0:
                PMR_wvf_S1[i]=np.random.normal(PMR_wvf_S1[i], sigmaAovA1*np.sqrt(PMR_wvf_S1[i]))
            if PMU_wvf_S1[i] > 0:
                PMU_wvf_S1[i]=np.random.normal(PMU_wvf_S1[i], sigmaAovA1*np.sqrt(PMU_wvf_S1[i]))
            if PMD_wvf_S1[i] > 0:
                PMD_wvf_S1[i]=np.random.normal(PMD_wvf_S1[i], sigmaAovA1*np.sqrt(PMD_wvf_S1[i]))
            
            
            if PML_wvf_S2[i] > 0:
                PML_wvf_S2[i]=np.random.normal(PML_wvf_S2[i], sigmaAovA1*np.sqrt(PML_wvf_S2[i]))
            if PMR_wvf_S2[i] > 0:
                PMR_wvf_S2[i]=np.random.normal(PMR_wvf_S2[i], sigmaAovA1*np.sqrt(PMR_wvf_S2[i]))
            if PMU_wvf_S2[i] > 0:
                PMU_wvf_S2[i]=np.random.normal(PMU_wvf_S2[i], sigmaAovA1*np.sqrt(PMU_wvf_S2[i]))
            if PMD_wvf_S2[i] > 0:
                PMD_wvf_S2[i]=np.random.normal(PMD_wvf_S2[i], sigmaAovA1*np.sqrt(PMD_wvf_S2[i]))
                
        
        #Generate waveforms (III. Convolute S2 with gap response -assumed to be a flat distribution)
        
        responseGap=np.ones(np.shape(Twvf))
        responseGap[Twvf>(Twvf[0]+gap/vdGap)]=0
        responseGap=1/Tbin*responseGap/(sum(responseGap))
        
        PML_wvf_S2=fftconvolve(PML_wvf_S2,responseGap, mode='full') #axes=Twvf)
        PMR_wvf_S2=fftconvolve(PMR_wvf_S2,responseGap, mode='full') #axes=Twvf)
        PMU_wvf_S2=fftconvolve(PMU_wvf_S2,responseGap, mode='full') #axes=Twvf)
        PMD_wvf_S2=fftconvolve(PMD_wvf_S2,responseGap, mode='full') #axes=Twvf)
        
        PML_wvf=PML_wvf_S1+PML_wvf_S2
        PMR_wvf=PMR_wvf_S1+PMR_wvf_S2
        PMU_wvf=PMU_wvf_S1+PMU_wvf_S2
        PMD_wvf=PMD_wvf_S1+PMD_wvf_S2
        
        #Generate waveforms (IV. Convolute with PM response)
        
        responsePM=np.exp(-(Twvf-meanTresponse)**2/(2*sigmaTresponse))
        responsePM=1/Tbin*responsePM/sum(responsePM)
        
        PML_wvf=fftconvolve(PML_wvf,responsePM, mode='full', axes=Twvf)
        PMR_wvf=fftconvolve(PMR_wvf,responsePM, mode='full', axes=Twvf)
        PMU_wvf=fftconvolve(PMU_wvf,responsePM, mode='full', axes=Twvf)
        PMD_wvf=fftconvolve(PMD_wvf,responsePM, mode='full', axes=Twvf)
        
        #Generate waveforms (V. Add noise, given in photon units)
        
        for i in range(len(Twvf)):
            
            PML_wvf[i]=PML_wvf[i]+np.random.noraml(0,sigmaNoise)
            PMR_wvf[i]=PMR_wvf[i]+np.random.noraml(0,sigmaNoise)
            PMU_wvf[i]=PMU_wvf[i]+np.random.noraml(0,sigmaNoise)
            PMD_wvf[i]=PMD_wvf[i]+np.random.noraml(0,sigmaNoise)
        
        #Generate waveforms (VI. Convolute with Osci response)
        
        Tau_osci=1/BW/2/pi
        responseOsci=1/Tau_osci*np.exp(-Twvf/Tau_osci)
        responseOsci=1/Tbin*responseOsci/sum(responseOsci)
        
        PML_wvf=fftconvolve(PML_wvf, responseOsci,mode='full', axes=twvf)
        PMR_wvf=fftconvolve(PMR_wvf, responseOsci,mode='full', axes=twvf)
        PMU_wvf=fftconvolve(PMU_wvf, responseOsci,mode='full', axes=twvf)
        PMD_wvf=fftconvolve(PMD_wvf, responseOsci,mode='full', axes=twvf)
        
        
        #Generate waveforms (VII. Convert from photon number to amplitude)
        
        PML_wvf=PML_wvf*A1
        PMR_wvf=PMR_wvf*A1
        PMU_wvf=PMU_wvf*A1
        PMD_wvf=PMD_wvf*A1
        
        
        return PML_wvf, PMR_wvf, PMU_wvf, PMD_wvf, Twvf, tPML_S1, tPMR_S1, tPMU_S1, tPMD_S1, tPML_S2, tPMR_S2, tPMU_S2, tPMD_S2







###     DEFINING THE INITIALIZATION PARAMS



OptFact = 10*10                  #Correction factor to adjust optical gain

P      = 0.1                     #ok<GPFST> %[bar] Operational conditions under CF4 in early 2019 (P.A. tfg). 
Epart  = 5.485                   # [MeV] Energy of alpha particles from 241Am (84.5% B.R. // 13% B.R., E = 5.442MeV).
Rpart  = 1.6/P * Epart / 5.4845  # [cm]  Range of alpha particles for CF4, according to TRIM (Manuel Caamaño). FIXME.
Mpart  = 3727                    # [MeV] 
z0     = 15                      # [cm]  The maximum drift distance of the TPC is 15cm.

nsteps = round(Rpart/0.005)      # Sample initial ionization/scintillation over 50um steps.
theta0 = 45                      # (deg). Default: 45
phi0   = 30                      # (deg). Default: 30

eps    = 1e-4                    # Just for numerics.










# IONIZATION


#global Wi, Ws, tau1, tau2, Ratio12, FanoQ, FanoS

Wi=54            #[ev] From "Properties of some gas mixtures used in tracking detectors", Archana Sharma, GSI Darmstadt, 1997
Ws=500           #[ev] At P~1 bar, and E<60V/cm/bar. From Morozov's: https://doi.org/10.1016/j.nima.2010.07.001 , https://doi.org/10.1016/j.nimb.2010.01.012.
tau1=2           #[ns] Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
tau2=40          #[ns] Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
Ratio12=0.35     #     Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
FanoQ=0.2        #     Typical Fano factor for ionization    (unknown)
FanoS=0.2        #     Typical Fano factor for scintillation (unknown)


# DRIFT


#global P, vd, DL, DT

vd = 3.5          # [cm/us]       Measured in Nausicaa1 at about 100 V/cm/bar (P. A. tfg), units changed below. 
DL = 250/np.sqrt(P)  # [um/sqrt(cm)] Simulated by Magboltz at about 100 V/cm/bar, units changed below.
DT = 250/np.sqrt(P)  # [um/sqrt(cm)] Simulated by Magboltz at about 100 V/cm/bar, units changed below.


# ANODE


#global hole, pitch, gap, OptGain, vdGap, tau1Gap, tau2Gap, Ratio12Gap

hole    = 0.05          # [cm] (diameter), from S. Williams / R. de Oliveira. 
pitch   = 0.1           # [cm] from S. Williams / R. de Oliveira (expected, but unconfirmed). 
gap     = 0.14          # [cm] (thickness of the structure).
OptGain = 105*OptFact   # Experimentally obtained with meshes in early 2019 (P.A. tfg).
vdGap   = 10            # [cm/us] preliminary simulated by Magboltz at about 10 kV/cm/bar (operating field in P.A. tfg was 25kV/cm/bar), units changed below.
tau1Gap = tau1          # Assume same spectrum for primary and secondary scintillation.
tau2Gap = tau2          # Assume same spectrum for primary and secondary scintillation.
Ratio12Gap = Ratio12    # Assume same spectrum for primary and secondary scintillation.


# PMT response 


#global zPM, rPM, RPM, QEPM, Tmesh, sigmaAovA1, sigmaNoise, A1, sigmaTresponse, meanTresponse, Tbin, BW

zPM            = 24        # [cm] PM z-position.
rPM            = 1.25      # [cm] PM radius.
RPM            = 7.2       # [cm] PM distance to center of flange.
QEPM           = 0.0885    # Effective quantum efficiency for PMs (P.A. tfg).
Tmesh          = 0.71      # Mesh optical transparency.
sigmaAovA1     = 1         # From Q analysis (M.F. tfg). Approximate.  
sigmaNoise     = 0.2       # In photon units (M.F. tfg). Approximate.
A1             = 5         # [mV] Single photon amplitude. Approximate.
sigmaTresponse = 10        # [ns] Width of PM time-response function. Approximate.
meanTresponse  = 20        # [ns] Delay for PM time-response function. Approximate.
Tbin           = 10        # [ns] Waveform sampling time.
BW             = 0.02      # [GHz] Scope bandwidth for waveform digitization.
if BW>0.1:                 # NOTE: it does not work for larger BW. This is not really a limitation, but it is strange.
    
    BW = 0.1  

# CMOS RESPONSE

#global M, lensN, sigmaNph, QECMOS, rebin, Npixel, pixelSize

M=-1/23              # magnification.
lensN=0.95           # lens number.
sigmaNph=1.2         # sigma in number of photons per pixel.
QECMOS=0.7*0.25      # 70% of the 25% of visible light.
rebin=10             # group pixels according to camera software.
Npixel=2048          # number of pixel per line.
pixelSize=6.5e-4     # pixel size [cm].








### ADJUST MAGNITUDES OF SWARM PARAMETERS

DL=DL*1e-4        #cm/sqrt(cm)
DT=DT*1e-4        #cm/sqrt(cm)
vd=vd*1e-3        #cm/ns
vdGap=vdGap*1e-3  #cm/ns




    