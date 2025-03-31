import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from numpy.linalg import eig

nrep =1 

f=open('input_data.dat','r')
lines=f.readlines()
for line in lines[0:]:
    data=line.split()
    if data[0]=='nsteps:':
        nsteps = int(data[1])
    if data[0]=='natoms:':
        natoms = int(data[1])
    if data[0]=='nonmvat:':
        nonmvat = int(data[1])
    if data[0]=='nmovingatoms:':
        nmovingatoms = int(data[1])
    if data[0]=='MDRestartFrequency:':
        MDRestartFrequency = int(data[1])
    if data[0]=='dt:':
        dt = float(data[1])
    if data[0]=='noh:':
        noh = int(data[1])
    if data[0]=='nwater:':
        nwater = int(data[1])
    if data[0]=='nox:':
        nox = int(data[1])
    if data[0]=='nhy:':
        nhy = int(data[1])
    if data[0]=='xbox:':
        xbox = float(data[1])
    if data[0]=='ybox:':
        ybox = float(data[1])
    if data[0]=='zbox:':
        zbox = float(data[1])

nsteps = int((nsteps/MDRestartFrequency)+1)-1 #+ int((nsteps2/MDRestartFrequency)+1)-2 #+int((nsteps3/MDRestartFrequency)+1)
#nsteps = 5000
#nsteps = 744
nstep = nsteps
        
print('nsteps:','  ',nsteps )
print('natoms:','  ',natoms)
print('nonmvat:','  ',nonmvat)
print('nmovingatoms:','  ',nmovingatoms)
print('MDRestartFrequency:','  ',MDRestartFrequency)
print('dt:','  ',dt)
print('noh:','  ',noh)
print('nwater:','  ',nwater)
print('nox:','  ',nox)
print('nhy:','  ',nhy)
print('xbox:','  ',xbox)
print('ybox:','  ',ybox)
print('zbox:','  ',zbox)

x_oh = np.empty(nsteps*nrep).reshape(nrep,nsteps)
y_oh = np.empty(nsteps*nrep).reshape(nrep,nsteps)
z_oh = np.empty(nsteps*nrep).reshape(nrep,nsteps)

kk = 0
while (kk<nrep):
    ii = 0
    f=open(f'oh_id_{ii+1}.dat','r')
    lines=f.readlines()
    for line in lines[:nsteps]:
        data=line.split()
        x_oh[kk,ii]=float(data[2])
        y_oh[kk,ii]=float(data[3])
        z_oh[kk,ii]=float(data[4])
        print(kk, ii, x_oh[kk,ii])
        ii = ii+1
    f.close()
    kk = kk+1

kk = 0
while (kk<nrep):
    ll = 0
    x_old = x_oh[kk,0]
    y_old = y_oh[kk,0]
    while ll<(nsteps):
        ijk = 1
        while ijk < 2:
            ijk = ijk+1
            if abs(x_oh[kk,ll]-x_old)>5:
                if x_oh[kk,ll]>x_old:
                    x_oh[kk,ll] = x_oh[kk,ll]-xbox
                    ijk = 1
                elif x_oh[kk,ll]<x_old:
                    x_oh[kk,ll] = x_oh[kk,ll]+xbox
                    ijk = 1
            if abs(y_oh[kk,ll]-y_old)>5: # atom need to be moved in y
                if y_oh[kk,ll]>y_old:
                    y_oh[kk,ll] = y_oh[kk,ll]-ybox
                    ijk = 1
                elif y_oh[kk,ll]<y_old:
                    y_oh[kk,ll] = y_oh[kk,ll]+ybox
                    ijk = 1
        x_old = x_oh[kk,ll]
        y_old = y_oh[kk,ll]
        ll = ll+1
    
    
    kk = kk+1

ndt = 500 #maximum tau value in fs/MDfreq
mmsdx = np.empty(ndt)
mmsdy = np.empty(ndt)
mmsdz = np.empty(ndt)
mmsdxy = np.empty(ndt)
mmsdxz = np.empty(ndt)
mmsdyz = np.empty(ndt)
dtime = np.empty(ndt)
dtime=  np.arange(0.01, (0.01*ndt+0.01), 0.01)
#MSD as a function of tau(time lag)


kk = 0
while kk<ndt:
    mmsdx[kk] = 0.0
    mmsdy[kk] = 0.0
    mmsdz[kk] = 0.0
    mmsdxy[kk] = 0.0
    mmsdxz[kk] = 0.0
    mmsdyz[kk] = 0.0
    jj = 0
    #print(((nsteps/(kk+1))-1))
    while jj< ((nsteps/(kk+1))-1):
        ll = 0
        while ll< (nrep):
           # print((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])**2)
            mmsdx[kk] = (x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])**2 +mmsdx[kk]
            mmsdy[kk] = (y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)])**2+mmsdy[kk]
            mmsdz[kk] = (z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)])**2+mmsdz[kk]
            mmsdxy[kk] = ((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])*(y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)]))+mmsdxy[kk]
            mmsdxz[kk] = ((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])*(z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)]))+mmsdxz[kk]
            mmsdyz[kk] = ((y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)])*(z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)]))+mmsdyz[kk]
            #print(kk,jj, mmsdx[kk])
            ll = ll+1
#        mmsdx[kk] = mmsdx[kk]/nrep
#        mmsdy[kk] = mmsdy[kk]/nrep
#        mmsdz[kk] = mmsdz[kk]/nrep
#        mmsdxy[kk] = mmsdxy[kk]/nrep
#        mmsdxz[kk] = mmsdxz[kk]/nrep
#        mmsdyz[kk] = mmsdyz[kk]/nrep
        jj = jj+1
 #   print(mmsdx[kk])
    mmsdx[kk] = mmsdx[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
#    print(kk, mmsdx[kk]) 
    mmsdy[kk] = mmsdy[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
    mmsdz[kk] = mmsdz[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
    mmsdxy[kk] = mmsdxy[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
    mmsdxz[kk] = mmsdxz[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
    mmsdyz[kk] = mmsdyz[kk]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
    kk = kk+1

fitx=np.polyfit(dtime[:],mmsdx[:],1)
fity=np.polyfit(dtime[:],mmsdy[:],1)
fitz=np.polyfit(dtime[:],mmsdz[:],1)
fitxy=np.polyfit(dtime[:],mmsdxy[:],1)
fitxz=np.polyfit(dtime[:],mmsdxz[:],1)
fityz=np.polyfit(dtime[:],mmsdyz[:],1)
#print('dtime and mmsdx')
#print(dtime)
#print(mmsdx)

pred_msdx=fitx[0]*dtime+fitx[1]
pred_msdy=fity[0]*dtime+fity[1]
pred_msdz=fitz[0]*dtime+fitz[1]

#calculating fitting error
errx = np.sqrt(np.sum((mmsdx-pred_msdx)**2)/len(mmsdx))
erry = np.sqrt(np.sum((mmsdy-pred_msdy)**2)/len(mmsdy))
errz = np.sqrt(np.sum((mmsdz-pred_msdz)**2)/len(mmsdz))
err2d = np.sqrt(errx**2+erry**2)
err2d
#find the optimum tau at lowest error

dxx = fitx[0]/2
dyy = fity[0]/2
dzz = fitz[0]/2
dxy = fitxy[0]/2
dxz = fitxz[0]/2
dyz = fityz[0]/2
#print(dxx)

dd2d = np.array([[dxx, dxy], 
                 [dxy, dyy]])
w2,v2=eig(dd2d)

print('Dx  =','{0: >#014.10f}'.format(w2[0]))
print('Dy  =','{0: >#014.10f}'.format(w2[1]))
print('D   =','{0: >#014.10f}'.format((w2[0]+w2[1])/2))
print('-----------------------')
print('Dx  =','{0: >#014.10f}'.format(dxx))
print('Dy  =','{0: >#014.10f}'.format(dyy))
print('D   =','{0: >#014.10f}'.format((dxx+dyy)/2))

#doh_data = open(f"doh_data_{nrep}.dat", "w")
#doh_data.write('nrep:' +str(nrep)+'\n')
#doh_data.write('RMSE at Tau= 5000 :'+str(err2d)+'\n')
#doh_data.write('Dx:'+ str(dxx)+'\n')
#doh_data.write('Dy:'+ str(dyy)+'\n')
#doh_data.write('D_tot:'+ str((dxx+dyy)/2))
#doh_data.close()
