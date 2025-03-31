import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from numpy.linalg import eig
import pandas as pd

def oh_xyz(nrep, nsteps):
    #nrep = number of replicas
    #nsteps = number of steps in fs/MD_freq (ex: for 10 ps simulation with MD_freq=10, nsteps=10,000/10=1000)

    #extracting x,y, z cordinates from oh_id.dat files
    x_oh = np.zeros((nrep,nsteps))
    y_oh = np.zeros((nrep,nsteps))
    z_oh = np.zeros((nrep,nsteps))

    for i in range(nrep):
        with open(f'oh_id_{i+1}.dat', 'r') as oh_id:
            xyz= oh_id.readlines()
            
            for j in range(len(xyz)):                
                x_oh[i,j]=float(xyz[j].split()[2])
                y_oh[i,j]=float(xyz[j].split()[3])
                z_oh[i,j]=float(xyz[j].split()[4])
                print(i,j, x_oh[i,j])
                
    #dtime=  np.arange(0.01, (0.01*ndt+0.01), 0.01)
    return(x_oh, y_oh, z_oh)

def oh_D(nrep, nsteps, natoms, nsheet, nwater, ncat, noh, xbox, ybox, zbox, ndt):
    #nrep = number of replicas
    #nsteps = number of steps in fs/MD_freq (ex: for 10 ps simulation with MD_freq=10, nsteps=10,000/10=1000)
    #natoms = total number of atoms in the sysytem
    #nsheet = total number of atoms in the graphene sheet
    #nwater = total number of water molecules
    #ncat   = total number of atoms in the cation/s
    #noh    =number of OH anions in the system
    #xbox, ybox, zbox = simulation box dimensions through x, y, z
    #ndt = maximum tau value in fs/MDfreq (ex: for Tau = 5000 fs with MD_freq=10, ndt=5000/10=500 )

    x_oh, y_oh, z_oh= oh_xyz(nrep, nsteps)
    
    MD_freq =10 #prints after every 10 steps (10 fs)
    dt =1.0 #in fs 
    
    nox = nwater + noh
    nhy = 2*nwater + noh
  

    msd_df = pd.DataFrame(columns=['Tau', 'dtime','msd_xx', 'msd_xy', 'msd_xz', 'msd_yy', 'msd_yz', 'msd_zz'])
    msd_df = msd_df.astype(float)
    

    dtime = np.empty(ndt)
    dtime=  np.arange(0.01, (0.01*ndt+0.01), 0.01)
    msdxx=np.zeros(ndt)
    msdxy=np.zeros(ndt)
    msdxz=np.zeros(ndt)
    msdyy=np.zeros(ndt)
    msdyz=np.zeros(ndt)
    msdzz=np.zeros(ndt)	

    
    for kk in range(ndt):
        msdxx = 0.0
        msdyy = 0.0
        msdzz = 0.0
        msdxy = 0.0
        msdxz = 0.0
        msdyz = 0.0
        
        for jj in range(math.ceil((nsteps/(kk+1))-1)): ##!!!recheck when kk=499
            
            for ll in range(nrep):
                msd_df.loc[kk, 'Tau']=kk+1
                msd_df.loc[kk, 'dtime']=(kk+1)/100
                #print((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])**2)
                msdxx = (x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])**2 +msdxx
                #print(kk, jj, msdxx)
                msdyy = (y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)])**2+msdyy
                msd_df.loc[kk, 'msd_yy']=msdyy
                msdzz = (z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)])**2+msdzz
                msd_df.loc[kk, 'msd_zz']=msdzz
                msdxy = ((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])*(y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)]))+msdxy
                msd_df.loc[kk, 'msd_xy']=msdxy
                msdxz = ((x_oh[ll,(jj*(kk+1)+kk+1)]-x_oh[ll,jj*(kk+1)])*(z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)]))+msdxz
                msd_df.loc[kk, 'msd_xz']=msdxz
                msdyz= ((y_oh[ll,(jj*(kk+1)+kk+1)]-y_oh[ll,jj*(kk+1)])*(z_oh[ll,(jj*(kk+1)+kk+1)]-z_oh[ll,jj*(kk+1)]))+msdyz
                msd_df.loc[kk, 'msd_yz']=msdyz
        
        msd_df.loc[kk, 'msd_xx']=msdxx
        #print(msd_df.loc[kk, 'msd_xx'])	
        # msd_df['msd_yx']=msd_df['msd_xy']
        # msd_df['msd_zx']=msd_df['msd_xz']
        # msd_df['msd_zy']=msd_df['msd_yz']
        msd_df.loc[kk, 'msd_xx'] = msd_df.loc[kk, 'msd_xx']/(((math.ceil(nsteps/(kk+1))-1))*nrep)
        #msd_df.iloc[kk, 2:] = msd_df.iloc[kk, 2:]/(((math.ceil(nsteps/(kk+1))-1))*nrep)
        #print(msd_df.loc[kk, 'msd_xx'])
       

    #perform linear fit
    fitxx=np.polyfit(msd_df['dtime'],msd_df['msd_xx'],1)
    fityy=np.polyfit(msd_df['dtime'],msd_df['msd_yy'],1)
    fitzz=np.polyfit(msd_df['dtime'],msd_df['msd_zz'],1)
    fitxy=np.polyfit(msd_df['dtime'],msd_df['msd_xy'],1)
    fitxz=np.polyfit(msd_df['dtime'],msd_df['msd_xz'],1)
    fityz=np.polyfit(msd_df['dtime'],msd_df['msd_yz'],1)
    #print('dtime and msd_xx')
    #print(msd_df['dtime'])
    #print(msd_df['msd_xx'])

    #create polynomial function
    pred_msd_xx = np.poly1d(fitxx)
    pred_msd_yy = np.poly1d(fityy)
    pred_msd_zz = np.poly1d(fitzz)

    #generate predicted values
    msd_df['pred_msd_xx'] = pred_msd_xx(msd_df['dtime'])
    msd_df['pred_msd_yy'] = pred_msd_yy(msd_df['dtime'])
    msd_df['pred_msd_zz'] = pred_msd_zz(msd_df['dtime'])

    #find the optimum tau at lowest error
    #calculate errors (Actual - Predicted)
    msd_df['xx_error'] = msd_df['msd_xx'] - msd_df['pred_msd_xx']  
    msd_df['yy_error'] = msd_df['msd_yy'] - msd_df['pred_msd_yy']  
    msd_df['zz_error'] = msd_df['msd_zz'] - msd_df['pred_msd_zz']  

    #MSE
    mse_xx = np.mean(msd_df['xx_error'] ** 2)
    mse_yy = np.mean(msd_df['yy_error'] ** 2)
    mse_zz = np.mean(msd_df['zz_error'] ** 2)
    
    #RMSE
    rmse_xx = np.sqrt(mse_xx)  
    rmse_yy = np.sqrt(mse_yy) 
    rmse_zz = np.sqrt(mse_zz) 
    rmse= np.sqrt(mse_xx+mse_yy) #considering x and y directions

    #generate 3D diffusion matrix
    dxx = fitxx[0]/2
    dyy = fityy[0]/2
    dzz = fitzz[0]/2
    dxy = fitxy[0]/2
    dxz = fitxz[0]/2
    dyz = fityz[0]/2
    print(dxx)
    d_3d = np.array([[dxx, dxy, dxz], 
                     [dxy, dyy, dyz],
                     [dxz, dyz, dzz]])
    w_3d,v_3d=eig(d_3d)

    #generate 2D diffusion matrix
    d_2d = np.array([[dxx, dxy], 
                     [dxy, dyy]])
    w_2d,v_2d=eig(d_2d)
    
    Dx  =dxx
    Dy  =dyy
    D_tot =(dxx+dyy)/2
    doh_data = open(f"doh_data_{nrep}.dat", "w")
    doh_data.write('nrep:' +str(nrep)+'\n')
    doh_data.write('RMSE at Tau= 5000 :'+str(rmse)+'\n')
    doh_data.write('Dx:'+ str(Dx)+'\n')
    doh_data.write('Dy:'+ str(Dy)+'\n')
    doh_data.write('D_tot:'+ str(D_tot))
    doh_data.close()

    #save msd_df
    msd_df.to_csv("msd_df.csv", index=False)
    
    return(rmse, Dx, Dy, D_tot)

#oh_D(nrep, nsteps, natoms, nsheet, nwater, ncat, noh, xbox, ybox, zbox, ndt)
rmse, Dx, Dy, D_tot= oh_D(1,      20000, 458,    287,    46,     31,   1,   15.198, 13.392, 15.0, 500)
print('RMSE at Tau= 5000 :',rmse)
print('Dx:', Dx)
print('Dy:', Dy)
print('D_total:', D_tot)


