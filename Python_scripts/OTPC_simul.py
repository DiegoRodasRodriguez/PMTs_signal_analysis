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

pi=math.pi


class OTPC_simul:
    
    
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
        
        #global: Wi, Ws, tau1, tau2, Ratio12, FanoQ, FanoS
        
        neToTmean=round(sum(DeltaE)/Wi)
        nphToTmean=round(sum(DeltaE)/Ws)
        
        sigmaStepQ=(np.sqrt(FanoQ)*np.sqrt(neToTmean))/np.sqrt(len(DeltaE))
        sigmaStepS=(np.sqrt(FanoS)*np.sqrt(nphToTmean))/np.sqrt(len(DeltaE))
        
        ne=np.zeros(len(DeltaE))
        nph=np.zeros(len(DeltaE))
        
        neToT=0
        nphToT=0
        
        for i in range(DeltaE):
            
            ne_mean=DeltaE[i]/Wi
            ne[i]=round(np.random.normal(ne_mean,sigmaStepQ))
            nph_mean=Delta[i]/Ws
            nph[i]=round(np.random.normal(nph_mean,sigmaStepS))
            neToT=neToT+ne[i]
            nphToT=nphToT+nph[i]
        
        xe=np.zeros(neToT)
        ye=np.zeros(neToT)
        ze=np.zeros(neToT)
        te=np.zeros(neToT)
        
        xph=np.zeros(nphToT)
        yph=np.zeros(nphToT)
        zph=np.zeros(nphToT)
        tph=np.zeros(nphToT)
        
        Lph=np.zeros(nphToT)
        
        I0=1
        nphCumm=0
        
        
        #Loop over energy deposits
        for i in range(len(DeltaE)):
            
            #electrons
            If=I0+ne[i]
            xe[I0:If]=x[i]
            ye[I0:If]=y[i]
            ze[I0:If]=z[i]
            te[I0:If]=t[i]
            I0=If+1
            
            #photons
            for j in range(nph[i]):
                
                rand_aux=np.random.rand()
                
                if rand_aux > Ratio12:
                    
                    tph[nphCumm+j]=np.random.exponential(tau2)
                    
                else:
                    tph[nphCumm+j]=np.random.exponential(tau1)
                
                
                xph[nphCumm+j]=x[i]
                yph[nphCumm+j]=y[i]
                zph[nphCumm+j]=z[i]
                Lph[nphCumm+j]=-1
                
            nphCumm=nphCumm+nph[i]
            
        return xe, ye, ze, te, xph, yph, zph, tph, Lph
    
    
    def TPCdrift(xe,ye,ze,te):
        
        #global: P, vd, DL, DT
        
        xeA=np.zeros(np.shape(xe))
        yeA=xeA
        zeA=xeA
        teA=xeA
        
        for i in range(len(xe)):
            
            teA[i]=te[i]+np.random.normal(ze[i]/vd,(DL*np.sqrt(ze[i]))/vd)
            xeA[i]=np.random.normal(xe[i],DT*np.sqrt(ze[i]))
            yeA=np.random.normal(ye[i],DT*np.sqrt(ze[i]))
            zeA[i]=0
        
        return xeA,yeA,zeA,teA
                    
    
    def TPCanode(xeA,yeA,zeA,teA):
        
        #global: hole, pitch, gap, OptGain, vdGap, tau1Gap, tau2Gap, Ratio12Gap
        
        xePh=pitch*round(xeA/pitch)
        yePh=pitch*round(yeA/pitch)
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
        
        #global: M, lensN, sigmaNph, QECMOS, rebin, Npixel, pixelSize Tmesh
        
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
        
        WePhMean=mean(WePh)     #Note: temporary patch
        
        for i in range(NpixXrebin):
            for j in range(NpixYrebin):
                
                NphXY[i,j]=np.random.poisson(WePhMean*WCMOS*NphXY[i,j])+np.random.normal(0,sigmaNph*np.sqrt(rebin*rebin))
        
        
        return NphXY
    
    
    
    def PMTresponse(xS2, yS2,zS2,tS2, WS2, xS1, yS1, zS1, tS1):


    
    