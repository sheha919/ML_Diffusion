import os
import numpy as np

def gen_therm_infiles(in_xyz, in_file, natoms, matoms, bx, by, bz):
    #This function is to generate input file for thermalization step using the given input geometry of the system (xyz file)

    #in_xyz     = initial geometry of the system (ex: "sys_geo" of sys_geo.xyz)
    #in_file    = name of the input file needs to generate (ex: dftb_in)
    #natoms     = number of total atoms
    #matoms     = number of moving atoms
    #bx, by, bz = dimensions of the box

    atom= []
    x=[]
    y=[]
    z= []
    with open(os.path.join(f"{in_xyz}.xyz"), 'r') as xyz_file:
        #for iqmol xyz
        xyz= xyz_file.readlines()[2:]
        #for regular xyz
    #    xyz= xyz_file.readlines()
        for indx,line in enumerate(xyz):
            atom.append(str(line.split()[0]))
            x.append(float(line.split()[1]))
            y.append(float(line.split()[2]))
            z.append(float(line.split()[3]))

    with open(os.path.join(f"{in_file}.hsd"), 'w') as in_file:
        in_file.write("#" + "\n" + "Geometry = GenFormat {" + "\n" + str(natoms) + " S" + "\n" + "C  H  O N" + "\n")

        for i in range(natoms):
            if atom[i]=='C':
                in_file.write(str(i+1) +"   "+ " 1" + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")
            if atom[i]=='H':
                in_file.write(str(i+1) +"   "+ " 2" + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")
            if atom[i]=='O':
                in_file.write(str(i+1) +"   "+ " 3" + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")
            if atom[i]=='N':
                in_file.write(str(i+1) +"   "+ " 4" + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")

        in_file.write("#dimensions of perodic box" + "\n")
        in_file.write( "     0.000     0.000     0.000" + "\n" + "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bz:
.3f}" + "\n" )
        in_file.write("}" + "\n" + "Driver = VelocityVerlet{" + "\n" + "  TimeStep [fs] = 1" + "\n")
        in_file.write("  Thermostat = NoseHoover {" + "\n" + "    Temperature [Kelvin] = 300" + "\n"+"    CouplingStrength [cm^-1] = 3700" + "\n"+ "  }"+ "\n")
        in_file.write("  Steps = 5000" + "\n")
        in_file.write("  MovedAtoms ="+ str(matoms) +":-1 # last -1 means all atoms; otherwise list first:last out of movable" + "\n" )
        in_file.write("  MDRestartFrequency = 10" + "\n" )
        in_file.write("  ConvergentForcesOnly = No #if the SCC cycle does not converge at any geometric step, forces will be calculated using the non-converged charges" + "\n" )
        in_file.write(" }" +"\n")
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


def gen_md_infiles(out_file, in_file, natoms, matoms, bx, by, bz):
    #This function is to extract the necessary information from output files and generate the next input file

    #out_file   = name of the output file produced after the simulation that used to extract the details for the next simulation (ex: geo_end)
    #in_file    = name of the input file needs to generate (ex: dftb_in)
    #natoms     = number of total atoms
    #matoms     = number of moving atoms
    #bx, by, bz = dimensions of the box

    atom= []
    x=[]
    y=[]
    z= []
    vx=[]
    vy=[]
    vz= []
    x_new=x.copy()
    y_new=y.copy()
    z_new= z.copy()
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

    #moving water molecules into the box
    for i in range(matoms+50,natoms+1,1):
        print(i)
        if abs(x[i])>bx/2:
            x_new[i]=(float(x[i]-bx))
        if abs(x[i])<bx/2:
            x_new[i]=(float(x[i]+bx))
        if abs(y[i])>by/2:
            y_new[i]=(float(y[i]-by))
        if abs(y[i])<by/2:
            y_new[i]=(float(y[i]+by))
        if abs(z[i])>bz/2:
            z_new[i]=(float(z[i]-bz))
        if abs(z[i])<bz/2:
            z_new[i]=(float(z[i]+bz))




    with open(os.path.join(f"{in_file}.hsd"), 'w') as in_file:
        in_file.write("#" + "\n" + "Geometry = GenFormat {" + "\n" + str(natoms) + " S" + "\n" + "C  H  O N" + "\n")

        for i in range(natoms):
            if atom[i]=='C':
                in_file.write(str(i+1) +"   "+ " 1" + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")
            if atom[i]=='H':
                in_file.write(str(i+1) +"   "+ " 2" + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")
            if atom[i]=='O':
                in_file.write(str(i+1) +"   "+ " 3" + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")
            if atom[i]=='N':
                in_file.write(str(i+1) +"   "+ " 4" + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")

        in_file.write("#dimensions of perodic box" + "\n")
        in_file.write( "     0.000     0.000     0.000" + "\n" + "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bz:
.3f}" + "\n" )
        in_file.write("}" + "\n" + "Driver = VelocityVerlet{" + "\n" + "  TimeStep [fs] = 1" + "\n" + "  Thermostat = None{}" + "\n" + "  Steps = 10000" + "\n")
        in_file.write("  MovedAtoms ="+ str(matoms) +":-1 # last -1 means all atoms; otherwise list first:last out of movable" + "\n" )
        in_file.write("  MDRestartFrequency = 10" + "\n" + "  ConvergentForcesOnly = Yes #if the SCC cycle does not converge at any geometric step, forces will be calculated using the non-converged charges" + "\
n" )
        in_file.write("  Velocities [AA/ps] {" + "\n" )

        for i in range(natoms):
            in_file.write("   " + f"{vx[i]:.8f}" + "   " + f"{vy[i]:.8f}" + "   " + f"{vz[i]:.8f}" +"\n")

        in_file.write(" }" +"\n")
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


#gen_therm_infiles('/home/shehani/DOE_project/int_xyz/c4_nmeth_double_cat_30h2o', 'dftb_in', 574, 431, 22.797, 13.392, 60)
#gen_md_infiles('geo_end', '../md1/dftb_in', 574, 431, 22.797, 13.392, 60)
gen_md_infiles('geo_end', 'test_in', 574, 431, 22.797, 13.392, 60)