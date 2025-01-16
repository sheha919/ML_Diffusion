#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np


# In[2]:


def gen_infiles(out_file, in_file, natoms, matoms, bx, by, bz):
    #This function is to extract the necessary information from output files and generate the next input file 

    atom= []
    x=[]
    y=[]
    z= []
    vx=[]
    vy=[]
    vz= []
    with open(os.path.join(f"{out_file}.xyz"), 'r') as out_file:
        xyz= out_file.readlines()
        end = len(xyz)
        start =end - natoms
        for i in range(start, end, 1):                
            atom.append(xyz[i].split()[0])
            x.append(float(xyz[i].split()[1]))
            y.append(float(xyz[i].split()[2]))
            z.append(float(xyz[i].split()[3]))
            vx.append(float(xyz[i].split()[5]))
            vy.append(float(xyz[i].split()[6]))
            vz.append(float(xyz[i].split()[7]))

    with open(os.path.join(f"{in_file}.hsd"), 'w') as in_file:
        in_file.write("#" + "\n" + "Geometry = GenFormat {" + "\n" + str(natoms) + " S" + "\n" + "C  H  O N" + "\n")

        for i in range(natoms):
            if atom[i]=='C':
                in_file.write(str(i+1) + " 1" + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i])+"\n")
            if atom[i]=='H':
                in_file.write(str(i+1) + " 2" + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i])+"\n")
            if atom[i]=='O':
                in_file.write(str(i+1) + " 3" + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i])+"\n")
            if atom[i]=='N':
                in_file.write(str(i+1) + " 4" + " " + str(x[i]) + " " + str(y[i]) + " " + str(z[i])+"\n")

        in_file.write("#dimensions of perodic box" + "\n")
        in_file.write( "     0.000     0.000     0.000" + "\n" + "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bz:.3f}" + "\n" )
        in_file.write("}" + "\n" + "Driver = VelocityVerlet{" + "\n" + "  TimeStep [fs] = 1" + "\n" + "  Thermostat = None{}" + "\n" + "  Steps = 10000" + "\n")  
        in_file.write("  MovedAtoms ="+ str(matoms) +":-1 # last -1 means all atoms; otherwise list first:last out of movable" + "\n" )
        in_file.write("  MDRestartFrequency = 10" + "\n" + "  ConvergentForcesOnly = Yes #if the SCC cycle does not converge at any geometric step, forces will be calculated using the non-converged charges" + "\n" )
        in_file.write("  Velocities [AA/ps] {" + "\n" )

        for i in range(natoms):
            in_file.write("   " + f"{vx[i]:.8f}" + "   " + f"{vy[i]:.8f}" + "   " + f"{vz[i]:.8f}" +"\n")

        in_file.write("	}" +"\n")             
        in_file.write("}" +"\n"+"\n" )           
        in_file.write("Hamiltonian = DFTB {"+"\n")
        in_file.write("  Scc = Yes" +"\n")
        in_file.write("  SCCTolerance = 1e-005"+"\n")
        in_file.write("  MaxSCCIterations= 500"+"\n")
        in_file.write("# set non-zero temperature for convergence"+"\n")
        in_file.write("  Filling = Fermi {"+"\n")
        in_file.write("  Temperature[K] = 300.0"+"\n")
        in_file.write("  }"+"\n"+"\n" )  
        in_file.write(" SlaterKosterFiles = Type2FileNames {"+"\n")
        in_file.write('  Prefix = "/home/shehani/DOE_project/slakos/" '+"\n")
        in_file.write('  Separator = "-"                     # Dash between type names'+"\n")
        in_file.write('  Suffix = ".skf"' +"\n")
        in_file.write(" }"+"\n"+"\n")
        in_file.write("  MaxAngularMomentum {"+"\n")
        in_file.write('    C = "p"'+"\n")
        in_file.write('    N = "p"'+"\n")
        in_file.write('    O = "p"'+"\n")
        in_file.write('    H = "s"'+"\n")
        in_file.write("  }"+"\n")
        in_file.write(" KPointsAndWeights = {"+"\n")
        in_file.write(" 0 0 0 1.0"+"\n")
        in_file.write(" }"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("Options {" +"\n")
        in_file.write(" WriteChargesAsText = yes"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("Analysis {"+"\n")
        in_file.write("  CalculateForces = Yes"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("ParserOptions {"+"\n")
        in_file.write("  ParserVersion = 7"+"\n")
        in_file.write("}"+"\n")
                    
         


# In[3]:


gen_infiles('geo_end', 'test', 335, 304, 15.198, 13.5, 40)


# In[ ]:




