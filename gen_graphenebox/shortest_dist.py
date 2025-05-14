import os
import numpy as np

def shortest_dist (file_name):
    #this function find the shortest distance between two atoms in two guest molecules
    #file_name: give the path to the xyz file of two guest molecules 
    file = os.path.join(f"{file_name}.xyz")
 
    indx=[]
    atom= []
    x=[]
    y=[]
    z= []
    with open(file, 'r') as file:
        xyz= file.readlines()
        nmol=int(len(xyz)/2)
        for i,line in enumerate(xyz):
            indx.append(int(i))
            atom.append(str(line.split()[0]))
            x.append(float(line.split()[1]))
            y.append(float(line.split()[2]))
            z.append(float(line.split()[3]))

    r =np.zeros((nmol, nmol))
    for i in range(nmol):
        for j in range(nmol, 2*nmol):
            rx = (x[i]-x[j])**2
            ry = (y[i]-y[j])**2
            rz = (z[i]-z[j])**2
            r[i,j-nmol] =np.sqrt(rx+ry+rz)
    short_r=np.min(r)
    ind = np.unravel_index(np.argmin(r), r.shape)
    indx1= indx[ind[0]]+1
    atom1=atom[ind[0]]
    indx2= indx[ind[1]+nmol]+1
    atom2=atom[ind[1]]

    print(f"shortest distance: between {atom1}({indx1}) and {atom2}({indx2}): {short_r} Angstroms  ")

    return(r)

#example
shortest_dist("BTD_A") 