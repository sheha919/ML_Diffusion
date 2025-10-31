import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd


def move_to_box(x, y, z, bx, by, bz):
    # insert x, y, z coordinates as 2D arrays and move atoms into the box
    a, b = x.shape

    x_new=x.copy()
    y_new=y.copy()
    z_new=z.copy()

    for i in range(a):
        for j in range(b):
            x_cor=x[i][j]
            if x_cor>bx/2:
                dx= int((abs(x_cor)-bx/2)/bx)+1
                x_new[i][j]=float(x_cor- (dx*bx))
            if x_cor<-bx/2:
                dx= int((abs(x_cor)-bx/2)/bx)+1
                x_new[i][j]=float(x_cor+ (dx*bx))

            y_cor=y[i][j]
            if (y_cor)>by/2:
                dy= int((abs(y_cor)-by/2)/by)+1
                y_new[i][j]=float(y_cor- (dy*by))
            if y_cor<-by/2:
                dy= int((abs(y_cor)-by/2)/by)+1
                y_new[i][j]=float(y_cor+ (dy*by))


            z_cor=z[i][j]
            if z_cor>bz/2:
                dz= int((abs(z_cor)-bz/2)/bz)+1
                z_new[i][j]=float(z_cor- (dz*bz))
            if z_cor<-bz/2:
                dz= int((abs(z_cor)-bz/2)/bz)+1
                z_new[i][j]=float(z_cor+ (dz*bz))

    return (x_new, y_new, z_new)

def dist_o_h(xo,yo,zo,xh,yh,zh):
    dist = np.sqrt((xh-xo)**2+(yh-yo)**2+(zh-zo)**2)
    return dist

def cal_msd(x1,y1,z1,x2,y2,z2):
    #calculates mean square displacement
    msd = np.array([0.0,0.0,0.0,0.0])
    msd[0] = ((x2-x1)**2 +(y2-y1)**2+(z2-z1)**2)
    msd[1] = (x2-x1)**2
    msd[2] = (y2-y1)**2
    msd[3] = (z2-z1)**2
    return msd

def gen_periodicbox(ox, x_bo, y_bo, z_bo, x_bh, y_bh, z_bh, xbox, ybox, zbox):
    #create 2D periodic box with only oxygen atoms in the system (creating 9 cells including the original box)
    #ox: indices of Oxygen atoms in the system
    #x_bo, y_bo, z_bo: arrays containing x, y, z coordinates after moving Oxygen atoms into the box
    #x_bh, y_bh, z_bh: arrays containing x, y, z coordinates after moving Hydrogen atoms into the box
    #xbox, ybox, zbox: Dimensions of the simulation box through x, y, z axes

    xperiodic_ox = np.append(x_bo,        (x_bo + (0*xbox)), axis=0) # original and 1
    xperiodic_ox = np.append(xperiodic_ox,(x_bo + (1*xbox)), axis=0) #2
    xperiodic_ox = np.append(xperiodic_ox,(x_bo + (1*xbox)),axis=0)  #3
    xperiodic_ox = np.append(xperiodic_ox,(x_bo + (1*xbox)),axis=0)  #4
    xperiodic_ox = np.append(xperiodic_ox,(x_bo + (0*xbox)),axis=0)  #5
    xperiodic_ox = np.append(xperiodic_ox,(x_bo +(-1*xbox)),axis=0)  #6
    xperiodic_ox = np.append(xperiodic_ox,(x_bo +(-1*xbox)),axis=0)  #7
    xperiodic_ox = np.append(xperiodic_ox,(x_bo +(-1*xbox)),axis=0)  #8

    yperiodic_ox = np.append(y_bo,        (y_bo + (1*ybox)),axis=0) # original and 1
    yperiodic_ox = np.append(yperiodic_ox,(y_bo + (1*ybox)),axis=0) #2
    yperiodic_ox = np.append(yperiodic_ox,(y_bo + (0*ybox)),axis=0) #3
    yperiodic_ox = np.append(yperiodic_ox,(y_bo +(-1*ybox)),axis=0) #4
    yperiodic_ox = np.append(yperiodic_ox,(y_bo +(-1*ybox)),axis=0) #5
    yperiodic_ox = np.append(yperiodic_ox,(y_bo +(-1*ybox)),axis=0) #6
    yperiodic_ox = np.append(yperiodic_ox,(y_bo + (0*ybox)),axis=0) #7
    yperiodic_ox = np.append(yperiodic_ox,(y_bo + (1*ybox)),axis=0) #8

    zperiodic_ox = (z_bo,z_bo,z_bo,z_bo,z_bo,z_bo,z_bo,z_bo,z_bo)
    zperiodic_ox = np.concatenate(zperiodic_ox)

    periodic_ox_indx=(ox,ox,ox,ox,ox,ox,ox,ox,ox)
    periodic_ox_indx=np.concatenate(periodic_ox_indx)

    return (periodic_ox_indx, xperiodic_ox, yperiodic_ox, zperiodic_ox)

def water(xo,yo,zo,xh,yh,zh,hy,periodic_ox_indx,nox,nhy):
    r_oh=np.zeros((nhy, nox*9))
    for i in range(r_oh.shape[0]): #go through Hs
        for j in range(r_oh.shape[1]):  #go through Os
            r_oh[i][j]=dist_o_h(xo[j],yo[j],zo[j],xh[i],yh[i],zh[i]) #distance between each O and H

    oh_indx = np.zeros((nhy,2), dtype=int) #for indices of O and H
    h2o = np.zeros((nox,6), dtype=int) #for water molecules

    for a in range(nhy):
        hh = r_oh[a][0] #distance with each H with first H
        oh_indx[a][0] = hy[a] #add index of each H
        oh_indx[a][1] = periodic_ox_indx[0] #add index of first Oxygen
        for b in range(1, nox*9): #go over in periodic box of Oxygen
            if r_oh[a][b]<hh: #find shotest bond
                oh_indx[a][1] = periodic_ox_indx[b]
                hh = r_oh[a][b]

    for nn in range(nox): #identify indices of atoms in water molecules
        h2o[nn][0]= periodic_ox_indx[nn] # add index of Oxygen to first column
        kk=1
        for mm in range(nhy):
            if oh_indx[mm][1]==periodic_ox_indx[nn]:
                h2o[nn][kk] = oh_indx[mm][0] #asign index of H to second/third column
                kk=kk+1
    return h2o #returns indices of O and H

def find_oh(h2o,nox):
    more_oh = 0

    with open('errors.dat','w') as err_file:

        for nn in range(nox):
            if h2o[nn][2]==0: #says it is OH-
                iioh = h2o[nn][0] #asign index of O
                more_oh = more_oh+1
            if more_oh>1:
                err_file.write('more than one oh')

            if h2o[nn][1]==0: #says it is O2-
                err_file.write('O atom found')

            if h2o[nn][3]>0: #says it is H3O+
                err_file.write('H3O+ found')

    return iioh #returns the index of O of OH-


def read_xyz(nsteps, natoms, nsheet, nwater, ncat, noh, xbox, ybox, zbox):
    #number of steps in fs/MD_freq (ex: for 10 ps simulation with MD_freq=10, nsteps=10,000/10=1000)
    #natoms = total number of atoms in the sysytem
    #nsheet = total number of atoms in the graphene sheet
    #nwater = total number of water molecules
    #ncat   = total number of atoms in the cation/s
    #noh    =number of OH anions in the system
    #xbox, ybox, zbox = simulation box dimensions through x, y, z

    MD_freq =10 #prints after every 10 steps (10 fs)
    dt =1.0 #in fs

    nox = nwater + noh
    nhy = 2*nwater + noh

    #assigning indices for moving Oxygen and Hydrogen atoms
    ox = []
    hy = []
    mv_start = nsheet + ncat +1
    mv_end   = natoms +1
    for i in range(mv_start, mv_end-1, 3):
        ox.append(i)
    for i in range(mv_start, mv_end, 1):
        if i not in ox:
            hy.append(i)

    #reading xyz coordinates related to Oxygens and Hydrogens
    x_h= np.zeros((nhy,nsteps))
    y_h= np.zeros((nhy,nsteps))
    z_h= np.zeros((nhy,nsteps))

    x_o= np.zeros((nox,nsteps))
    y_o= np.zeros((nox,nsteps))
    z_o= np.zeros((nox,nsteps))


    with open('geo_end.xyz', 'r') as out_file:
        xyz= out_file.readlines()
        for i in range(nsteps):
            end = (i+1)*(natoms+2)
            start= end - natoms
            io=0
            ih=0
            for j in ox:
                oindx = start+j-1
                x_o[io][i]=float(xyz[oindx].split()[1])
                y_o[io][i]=float(xyz[oindx].split()[2])
                z_o[io][i]=float(xyz[oindx].split()[3])
                io=io+1
            for k in hy:
                hindx = start+k-1
                x_h[ih][i]=float(xyz[hindx].split()[1])
                y_h[ih][i]=float(xyz[hindx].split()[2])
                z_h[ih][i]=float(xyz[hindx].split()[3])
                ih = ih+1
    print('done reading coordinates')
    return (x_o, y_o, z_o,x_h, y_h, z_h, ox, hy)

def oh_identifier(nsteps, natoms, nsheet, nwater, ncat, noh, xbox, ybox, zbox):
    #number of steps in fs/MD_freq (ex: for 10 ps simulation with MD_freq=10, nsteps=10,000/10=1000)
    #natoms = total number of atoms in the sysytem
    #nsheet = total number of atoms in the graphene sheet
    #nwater = total number of water molecules
    #ncat   = total number of atoms in the cation/s
    #noh    =number of OH anions in the system
    #xbox, ybox, zbox = simulation box dimensions through x, y, z

    MD_freq =10 #prints after every 10 steps (10 fs)
    dt =1.0 #in fs

    nox = nwater + noh
    nhy = 2*nwater + noh

    #assigning indices for moving Oxygen and Hydrogen atoms
    ox = []
    hy = []
    mv_start = nsheet + ncat +1
    mv_end   = natoms +1
    for i in range(mv_start, mv_end-1, 3):
        ox.append(i)
    for i in range(mv_start, mv_end, 1):
        if i not in ox:
            hy.append(i)

    #reading xyz coordinates related to Oxygens and Hydrogens
    x_h= np.zeros((nhy,nsteps))
    y_h= np.zeros((nhy,nsteps))
    z_h= np.zeros((nhy,nsteps))

    x_o= np.zeros((nox,nsteps))
    y_o= np.zeros((nox,nsteps))
    z_o= np.zeros((nox,nsteps))


    with open('geo_end.xyz', 'r') as out_file:
        xyz= out_file.readlines()
        for i in range(nsteps):
            end = (i+1)*(natoms+2)
            start= end - natoms
            io=0
            ih=0
            for j in ox:
                oindx = start+j-1
                x_o[io][i]=float(xyz[oindx].split()[1])
                y_o[io][i]=float(xyz[oindx].split()[2])
                z_o[io][i]=float(xyz[oindx].split()[3])
                io=io+1
            for k in hy:
                hindx = start+k-1
                x_h[ih][i]=float(xyz[hindx].split()[1])
                y_h[ih][i]=float(xyz[hindx].split()[2])
                z_h[ih][i]=float(xyz[hindx].split()[3])
                ih = ih+1
    print('done reading coordinates')

    x_bo, y_bo, z_bo = move_to_box(x_o, y_o, z_o, xbox, ybox, zbox)
    x_bh, y_bh, z_bh = move_to_box(x_h, y_h, z_h, xbox, ybox, zbox)
    print('done moving coordinates to the box')

    periodic_ox_indx, xperiodic_ox, yperiodic_ox, zperiodic_ox = gen_periodicbox(ox, x_bo, y_bo, z_bo, x_bh, y_bh, z_bh,xbox, ybox, zbox)
    print('done making the periodic cells')

    x_oh = np.zeros(nsteps)
    y_oh = np.zeros(nsteps)
    z_oh = np.zeros(nsteps)
    indx_oh = np.zeros(nsteps)

    msd_oh=np.zeros(nsteps)
    msd_ohx=np.zeros(nsteps)
    msd_ohy=np.zeros(nsteps)
    msd_ohz=np.zeros(nsteps)

    xyz_oh = pd.DataFrame()
    rattle = open("rattle.dat", "w")
    for ll in range(nsteps):
        xo = [sub[ll] for sub in xperiodic_ox]
        yo = [sub[ll] for sub in yperiodic_ox]
        zo = [sub[ll] for sub in zperiodic_ox]
        xh = [sub[ll] for sub in x_bh]
        yh = [sub[ll] for sub in y_bh]
        zh = [sub[ll] for sub in z_bh]
        h2o = water(xo,yo,zo,xh,yh,zh,hy,periodic_ox_indx,nox,nhy)
        indx_o = find_oh(h2o,nox) #index of O of OH-
        indx_ox =ox.index(indx_o) #index of ox array corresponding to real indx_o
        if ll==0:
            xo_0 = x_o[indx_ox][0]
            yo_0 = y_o[indx_ox][0]
            zo_0 = z_o[indx_ox][0]
            xo_int= xo_0
            yo_int= yo_0

        xmv = 'no change'
        ymv = 'no change'

        #checking periodic conditions
        # atom needs to be moved to the box if it moves 5 Angstroms within 10 fs relative to the previous time step
        distx= abs(x_o[indx_ox][ll]-xo_int)
        disty= abs(y_o[indx_ox][ll]-yo_int)

        if distx>5: # atom need to be moved in 5 Angstroms within 10 fs
            print(ll, 'distx')
            dxx=(distx/xbox)
            dxxi=int(distx/xbox)
            if (dxx-dxxi)>0.8:
                dxxi=dxxi+1
            if x_o[indx_ox][ll]>xo_int:
                x_o[indx_ox][ll]= x_o[indx_ox][ll]-(xbox*dxxi)
                xmv = 'atom moved in negative x'
            elif x_o[indx_ox][ll]<xo_int:
                x_o[indx_ox][ll]= x_o[indx_ox][ll]+(xbox*dxxi)
                xmv = 'atom moved in positive x'

        if disty>5: # atom need to be moved in y
            print(ll, 'disty')
            dyy=(disty/ybox)
            dyyi=int(disty/ybox)
            if (dyy-dyyi)>0.8:
                 dyyi=dyyi+1
            if y_o[indx_ox][ll]>yo_int:
                y_o[indx_ox][ll]= y_o[indx_ox][ll]-(ybox*dyyi)
                ymv = 'atom moved in negative y'
            elif y_o[indx_ox][ll]<yo_int:
                y_o[indx_ox][ll]= y_o[indx_ox][ll]+(ybox*dyyi)
                ymv = 'atom moved in positive y'


        xo_int = x_o[indx_ox][ll]
        yo_int = y_o[indx_ox][ll]

        #extract x, y, z coordinates of O of OH- in each step
        x_oh[ll] = x_o[indx_ox][ll]
        y_oh[ll] = y_o[indx_ox][ll]
        z_oh[ll] = z_o[indx_ox][ll]

        #add index of O of OH- in each step
        indx_oh[ll] = indx_o


        #remove rattling
        if (indx_oh[ll]==indx_oh[(ll-20)]) and (indx_oh[ll]!=indx_oh[(ll-1)]): #if rattling occurs for 200 fs
            xmv = 'no change'
            ymv = 'no change'
            indx_rox =ox.index(indx_oh[ll]) #index of ox array corresponding to real indx_oh[ll]
            rattle.write('rattling1 starts here'+ '\n')
            print(ll, 'rat1')
            for llk in range((ll-19),ll):
                x_oh[llk] = x_o[indx_rox][llk]
                y_oh[llk] = y_o[indx_rox][llk]
                z_oh[llk] = z_o[indx_rox][llk]
                indx_oh[llk] = indx_oh[ll]


                distrx=abs(x_oh[llk]-x_oh[llk-1])
                distry=abs(y_oh[llk]-y_oh[llk-1])


                if distrx>5: # atom need to be moved in x
                    drxx=(distrx/xbox)
                    drxxi=int(distrx/xbox)
                    print(ll, llk, 'distrx')
                    if (drxx-drxxi)>0.8:
                        drxxi=drxxi+1
                    if x_oh[llk]>x_oh[llk-1]:
                        x_oh[llk] = x_oh[llk]-(xbox*drxxi)
                        xmv = 'atom moved in negative x'
                    elif x_oh[llk]<x_oh[llk-1]:
                        x_oh[llk] = x_oh[llk]+(xbox*drxxi)
                        xmv = 'atom moved in positive x'

                if distry>5: # atom need to be moved in y
                    dryy=(distry/ybox)
                    dryyi=int(distry/ybox)
                    print(ll, llk, 'distry')
                    if (dryy-dryyi)>0.8:
                        dryyi= dryyi+1
                    if y_oh[llk]>y_oh[llk-1]:
                        y_oh[llk] = y_oh[llk]-(ybox*dryyi)
                        ymv = 'atom moved in negative y'
                    elif y_oh[llk]<y_oh[llk-1]:
                        y_oh[llk] = y_oh[llk]+(ybox*dryyi)
                        ymv = 'atom moved in positive y'

                rattle.write('{:<6d}'.format(llk)+'{:<6d}'.format(int(indx_oh[llk]))+'{0: >#016.8f}'.format(x_oh[llk])+'{0: >#016.8f}'.format(y_oh[llk])+'{0: >#016.8f}'.format(z_oh[llk])+'  '+xmv+'  '+ymv+ '\n')


            xo_int = x_oh[llk+1]
            yo_int = y_oh[llk+1]
            rattle.write('rattling1 ends here'+ '\n')


        # remove rattling just after hopping (50 fs)
        if (indx_oh[ll]==indx_oh[(ll-5)]) and (indx_oh[ll]!=indx_oh[(ll-1)]):
            xmv = 'no change'
            ymv = 'no change'
            indx_hox =ox.index(indx_oh[ll]) #index of ox array corresponding to real indx_oh[ll]
            rattle.write('rattling2 starts here'+ '\n')
            print(ll, 'rat2')

            for llk in range((ll-4), ll):
                x_oh[llk] = x_o[indx_hox][llk]
                y_oh[llk] = y_o[indx_hox][llk]
                z_oh[llk] = z_o[indx_hox][llk]
                indx_oh[llk] = indx_oh[ll]


                disthx=abs(x_oh[llk]-x_oh[llk-1])
                disthy=abs(y_oh[llk]-y_oh[llk-1])

                if disthx>5: #atoms need to be moved to the box through x
                    dhxx=(disthx/xbox)
                    dhxxi=int(disthx/xbox)
                    print(ll, llk, 'disthx')
                    if (dhxx-dhxxi)>0.8:
                        dhxxi=dhxxi+1
                    if x_oh[llk]>x_oh[llk-1]:
                        x_oh[llk] = x_oh[llk]-(xbox*dhxxi)
                        xmv = 'atom moved in negative x'
                    elif x_oh[llk]<x_oh[llk-1]:
                        x_oh[llk] = x_oh[llk]+(xbox*dhxxi)
                        xmv = 'atom moved in positive x'

                if disthy>5: #atoms need to be moved to the box through y
                    dhyy=(disthy/ybox)
                    dhyyi=int(disthy/ybox)
                    print(ll, llk, 'disthy')
                    if (dhyy-dhyyi)>0.8:
                        dhyyi=dhyyi+1
                    if y_oh[llk]>y_oh[llk-1]:
                        y_oh[llk] = y_oh[llk]-(ybox*dhyyi)
                        ymv = 'atom moved in negative y'
                    elif y_oh[llk]<y_oh[llk-1]:
                        y_oh[llk] = y_oh[llk]+(ybox*dhyyi)
                        xyz_oh.loc[llk, "y"] =xyz_oh.loc[llk, "y"]+(ybox*dhyyi)
                        ymv = 'atom moved in positive y'

                rattle.write('{:<6d}'.format(llk)+'{:<6d}'.format(int(indx_oh[llk]))+'{0: >#016.8f}'.format(x_oh[llk])+'{0: >#016.8f}'.format(y_oh[llk])+'{0: >#016.8f}'.format(z_oh[llk])+'  '+xmv+'  '+ymv+ '\n')

            xo_int = x_oh[llk+1]
            yo_int = y_oh[llk+1]
            rattle.write('rattling2 ends here'+ '\n')

        distx_n= abs(x_oh[ll] - x_oh[ll-1])
        distx_npos=abs(x_oh[ll] + xbox - x_oh[ll-1])
        distx_nneg=abs(x_oh[ll] - xbox - x_oh[ll-1])

        disty_n= abs(y_oh[ll]- y_oh[ll-1])
        disty_npos= abs(y_oh[ll] + ybox - y_oh[ll-1])
        disty_nneg= abs(y_oh[ll] - ybox - y_oh[ll-1])

        x_dist= min(distx_n, distx_npos, distx_nneg)
        if x_dist== distx_n:
            x_oh[ll]=x_oh[ll]
        elif x_dist== distx_npos:
            x_oh[ll]=x_oh[ll]+xbox
        elif x_dist== distx_nneg:
            x_oh[ll]=x_oh[ll]-xbox
        xo_int=x_oh[ll]

        y_dist= min(disty_n, disty_npos, disty_nneg)
        if y_dist== disty_n:
            y_oh[ll]=y_oh[ll]
        elif y_dist== disty_npos:
            y_oh[ll]=y_oh[ll]+ybox
        elif y_dist== disty_nneg:
            y_oh[ll]=y_oh[ll]-ybox
        yo_int=y_oh[ll]


    rattle.close()
    print("done oh_identifier")

    return(xo_0,yo_0,zo_0,indx_oh, x_oh,y_oh,z_oh)


#######################################################################################################################################################################################################
#(nsteps, natoms, nsheet, nwater, ncat, noh, xbox, ybox, zbox)
oh_x0,oh_y0, oh_z0, indx, x, y, z=oh_identifier(5001, 617, 431, 55, 19, 1, 22.797, 13.392, 12.0)
#xo_01, yo_01, zo_01, indx_oh1, x_oh1,y_oh1,z_oh1



#md id
n_md=1

oh_msd = pd.DataFrame(columns=['nstep', 'msd_oh', 'msd_ohx', 'msd_ohy', 'msd_ohz'])
nsteps=len(x)

with open(f'oh.{n_md}_id.dat', 'w') as oh_id:
    for i in range(nsteps):

        oh_id.write('{:<6d}'.format(i) +'{:<6d}'.format(int(indx[i]))+'{0: >#016.8f}'.format(x[i])+'{0: >#016.8f}'.format(y[i])+'{0: >#016.8f}'.format(z[i])+'\n')
#oh_msd.to_csv(path+"oh_msd.csv", index=False)

plt.figure()
plt.plot(x[:],color='black',linewidth=2,label="X")
plt.plot(y[:],color='r',linewidth=2,label="Y")
plt.plot(z[:],color='g',linewidth=2,label="Z")
plt.legend(loc="upper left")
plt.xlabel('Steps')
plt.ylabel('position')
plt.title("OH Diffusion")
plt.tick_params(axis="x",which='major', direction="in", length=5, width=0.5)
plt.tick_params(axis="y",which='major',direction="in", length=5, width=0.5)
plt.tick_params(axis="x",which='minor', direction="in", length=5, width=0.5)
plt.tick_params(axis="y",which='minor',direction="in", length=5, width=0.5)
plt.grid(which='major',color='#CCCCCC', linestyle='--', linewidth=0.5)
plt.ticklabel_format(axis="x", style="plain", scilimits=(0,0),useMathText=True)
#plt.ylim([0,80])
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.savefig(f'oh.{n_md}.jpg', dpi=100)