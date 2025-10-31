import os
import numpy as np


def gen_therm_infiles(in_xyz, in_file, natoms, matoms, bx, by, bz):
    #This function is to generate input file for thermalization step using the given input geometry of the system (xyz file)

    #in_xyz     = initial geometry of the system (ex: "sys_geo" of sys_geo.xyz)
    #in_file    = name of the input file needs to generate (ex: dftb_in)
    #natoms     = number of total atoms
    #matoms     = starting index of moving atoms
    #bx, by, bz = dimensions of the box

    #alowed not to interact through z-axis
    bZ= 45+ bz

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
        in_file.write( "     0.000     0.000     0.000"+ "\n")
        in_file.write( "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bZ:.3f}" + "\n" )
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
        in_file.write('  Prefix = "/home/shehani/DOE_project/slakos_3ob_nh/" '+"\n")
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
        in_file.write(" TimingVerbosity = 2"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("Analysis {"+"\n")
        in_file.write("  CalculateForces = Yes"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("ParserOptions {"+"\n")
        in_file.write("  ParserVersion = 7"+"\n")
        in_file.write("}"+"\n")


def gen_md_infiles(out_file, in_file, xyz_path, natoms, matoms, bx, by, bz):
    #This function is to extract the necessary information from output files and generate the next input file

    #out_file   = name of the output file produced after the simulation that used to extract the details for the next simulation (ex: geo_end)
    #in_file    = name of the input file needs to generate (ex: dftb_in)
    #natoms     = number of total atoms
    #matoms     = starting index of moving atoms
    #cat_atoms  = number of cations in single cation
    #ncat       = number of cations
    #bx, by, bz = dimensions of the box

    #alowed not to interact through z-axis
    bZ= 45+ bz

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

    #generating xyz file with last step
    with open(os.path.join(f"{xyz_path}/final_geo.xyz"), 'w') as final_xyz:
        final_xyz.write(str(natoms) + "\n" + "\n")
        for i in range(natoms):
            final_xyz.write(atom[i] + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")

    x_new=x.copy()
    y_new=y.copy()
    z_new= z.copy()

    #moving water molecules into the box
    for i in range(matoms-1, natoms,1):
        if (x[i])>bx/2:
            dx= int((abs(x[i])-bx/2)/bx)+1
            x_new[i]=float(x[i]-(dx*bx))
        if (x[i])<-bx/2:
            dx= int((abs(x[i])-bx/2)/bx)+1
            x_new[i]=float(x[i]+ (dx*bx))
        if (y[i])>by/2:
            dy= int((abs(y[i])-by/2)/by)+1
            y_new[i]=float(y[i]- (dy*by))
        if (y[i])<-by/2:
            dy= int((abs(y[i])-by/2)/by)+1
            y_new[i]=float(y[i]+ (dy*by))
        if (z[i])>bz/2:
            dz= int((abs(z[i])-bz/2)/bz)+1
            z_new[i]=float(z[i]- (dz*bz))
        if (z[i])<-bz/2:
            dz= int((abs(z[i])-bz/2)/bz)+1
            z_new[i]=float(z[i]+ (dz*bz))

    #generating xyz file after moving water and OH into the box
    with open(os.path.join(f"{xyz_path}/movedin.xyz"), 'w') as mv_file:
        mv_file.write(str(natoms) + "\n" + "\n")
        for i in range(natoms):
            mv_file.write(atom[i] + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")

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
        in_file.write( "     0.000     0.000     0.000"+ "\n")
        in_file.write( "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bZ:.3f}" + "\n" )
        in_file.write("}" + "\n" + "Driver = VelocityVerlet{" + "\n" + "  TimeStep [fs] = 1" + "\n" + "  Thermostat = None{}" + "\n" + "  Steps = 200000" + "\n")
        in_file.write("  MovedAtoms ="+ str(matoms) +":-1 # last -1 means all atoms; otherwise list first:last out of movable" + "\n" )
        in_file.write("  MDRestartFrequency = 10" + "\n" + "  ConvergentForcesOnly = Yes #if the SCC cycle does not converge at any geometric step, forces will be calculated using the non-converged charges" + "\n" )
        in_file.write("  Velocities [AA/ps] {" + "\n" )

        for i in range(natoms):
            in_file.write("   " + f"{vx[i]:.8f}" + "   " + f"{vy[i]:.8f}" + "   " + f"{vz[i]:.8f}" +"\n")

        in_file.write(" }" +"\n")
        in_file.write("}" +"\n"+"\n" )
        in_file.write("Hamiltonian = DFTB {"+"\n")
        in_file.write("  Scc = Yes" +"\n")
        in_file.write("  SCCTolerance = 1e-001"+"\n")
        in_file.write("  MaxSCCIterations= 500"+"\n")
        in_file.write("# set non-zero temperature for convergence"+"\n")
        in_file.write("  Filling = Fermi {"+"\n")
        in_file.write("  Temperature[K] = 300.0"+"\n")
        in_file.write("  }"+"\n"+"\n" )
        in_file.write(" SlaterKosterFiles = Type2FileNames {"+"\n")
        in_file.write('  Prefix = "/home/shehani/DOE_project/slakos_3ob_nh/" '+"\n")
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
        in_file.write(" TimingVerbosity = 2"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("Analysis {"+"\n")
        in_file.write("  CalculateForces = Yes"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("ParserOptions {"+"\n")
        in_file.write("  ParserVersion = 7"+"\n")
        in_file.write("}"+"\n")


def gen_md_replica(cat_sys, out_file, in_file, xyz_path, istep, natoms, matoms, bx, by, bz):
    #This function is to extract the necessary information from output files and generate the input file for replicas
    #cat_sys = file containing coordinates for two graphene sheets with cation. This is same for all replicas in a given cation
    #istep is the cycle of the MD step i.e (istep 0==>MD step 0, itep 1==> MD step 10)
    #out_file   = name of the output file produced after the simulation that used to extract the details for the next simulation (ex: geo_end)
    #in_file    = name of the input file needs to generate (ex: dftb_in)
    #xyz_path   = path to extracted xyz files
    #natoms     = number of total atoms
    #matoms     = starting index of moving atoms
    #bx, by, bz = dimensions of the box

    #alowed not to interact through z-axis
    bZ= 45+ bz

    atom= []
    x=[]
    y=[]
    z= []
   #extract coordinates for graphene sheets with cation==> does not change for all replicas
    with open(os.path.join(f"/home/shehani/DOE_project/int_xyz/{cat_sys}.xyz"), 'r') as sys:
        sys_xyz= sys.readlines()[2:]
        for indx,line in enumerate(sys_xyz):
            atom.append(str(line.split()[0]))
            x.append(float(line.split()[1]))
            y.append(float(line.split()[2]))
            z.append(float(line.split()[3]))

    with open(os.path.join(f"{out_file}.xyz"), 'r') as out_file:
        xyz= out_file.readlines()[2:]
        for indx, line in enumerate(xyz):
            start = istep*(natoms+2)+ len(sys_xyz)
            end = (istep*(natoms+2)) + natoms

        for i in range(start, end, 1):
            atom.append(xyz[i].split()[0])
            x.append(float(xyz[i].split()[1]))
            y.append(float(xyz[i].split()[2]))
            z.append(float(xyz[i].split()[3]))
    #generating xyz file with istep
    with open(os.path.join(f"{xyz_path}/{istep}_step_geo.xyz"), 'w') as ext_xyz:
        ext_xyz.write(str(natoms) + "\n" + "\n")
        for i in range(natoms):
            ext_xyz.write(atom[i] + "   " + f"{x[i]:.8f}" + "   " + f"{y[i]:.8f}" + "   " + f"{z[i]:.8f}"+"\n")

    x_new=x.copy()
    y_new=y.copy()
    z_new= z.copy()

    #moving water molecules into the box
    for i in range(matoms-1, natoms,1):
        if (x[i])>bx/2:
            dx= int((abs(x[i])-bx/2)/bx)+1
            x_new[i]=float(x[i]-(dx*bx))
        if (x[i])<-bx/2:
            dx= int((abs(x[i])-bx/2)/bx)+1
            x_new[i]=float(x[i]+ (dx*bx))
        if (y[i])>by/2:
            dy= int((abs(y[i])-by/2)/by)+1
            y_new[i]=float(y[i]- (dy*by))
        if (y[i])<-by/2:
            dy= int((abs(y[i])-by/2)/by)+1
            y_new[i]=float(y[i]+ (dy*by))
        if (z[i])>bz/2:
            dz= int((abs(z[i])-bz/2)/bz)+1
            z_new[i]=float(z[i]- (dz*bz))
        if (z[i])<-bz/2:
            dz= int((abs(z[i])-bz/2)/bz)+1
            z_new[i]=float(z[i]+ (dz*bz))

    #generating xyz file after moving water and OH into the box
    with open(os.path.join(f"{xyz_path}/{istep}_step_moved.xyz"), 'w') as mv_file:
        mv_file.write(str(natoms) + "\n" + "\n")
        for i in range(natoms):
            mv_file.write(atom[i] + "   " + f"{x_new[i]:.8f}" + "   " + f"{y_new[i]:.8f}" + "   " + f"{z_new[i]:.8f}"+"\n")

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
        in_file.write( "     0.000     0.000     0.000"+ "\n")
        in_file.write( "    " + f"{bx:.3f}" + "     0.000     0.000" + "\n" + "     0.000     " + f"{by:.3f}" + "    0.000" + "\n" + "     0.000     0.000     " + f"{bZ:.3f}" + "\n" )
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
        in_file.write('  Prefix = "/home/shehani/DOE_project/slakos_3ob_nh/" '+"\n")
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
        in_file.write(" TimingVerbosity = 2"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("Analysis {"+"\n")
        in_file.write("  CalculateForces = Yes"+"\n")
        in_file.write("}"+"\n"+"\n")
        in_file.write("ParserOptions {"+"\n")
        in_file.write("  ParserVersion = 7"+"\n")
        in_file.write("}"+"\n")

####################################################################change start and end############################################################################################
st=1
end=5
for i in range(st,end+1,1):
        if not os.path.exists(f"therm{i}"):

        #gen_md_infiles(f"therm{i}/geo_end", f"md{i}/dftb_in", f"therm{i}", 449, 288, 15.198, 13.392, 15)
####################################################################change parameters##################################################################################################
                os.makedirs(f"therm{i}")
        gen_md_replica('c4c4_dualcat_z18','/home/shehani/DOE_project/frozen_cat/double_cat/c4c4_98h2o_inv/therm1/geo_end',f"therm{i}/dftb_in",f"therm{i}",i*100, 922, 575, 30.396, 13.392, 18)

        with open(os.path.join(f"therm{i}/sub_dftb.sh"), 'w') as sub:
                sub.write("#!/bin/sh" +  "\n")
####################################################################change job name##################################################################################################
                sub.write(f"#SBATCH --job-name=c4c4.{i}_inv" + "\n")
                with open(os.path.join("sub_dftb.sh"), 'r') as or_sub:
                        lines = or_sub.readlines()[2:]
                        for line in lines:
                                sub.write(line)

        with open(os.path.join(f"therm{i}/therm_analysis.py"), 'w') as anal:
                with open(os.path.join("therm_analysis.py"), 'r') as or_anal:
                        lines = or_anal.readlines()
                        for line in lines:
                                anal.write(line)
####################################################################change numer of water and steps ##################################################################################################
                anal.write("md_analysis(98, 4000, 5000)")


        with open(os.path.join(f"therm{i}/python.sh"), 'w') as py1:
                with open(os.path.join("python.sh"), 'r') as or_py1:
                        lines = or_py1.readlines()
                        for line in lines:
                                py1.write(line)

        with open(os.path.join(f"therm{i}/py_therm.sh"), 'w') as py2:
                with open(os.path.join("py_therm.sh"), 'r') as or_py2:
                        lines = or_py2.readlines()
                        for line in lines:
                                py2.write(line)

        with open(os.path.join("sub_multi.sh"), 'w') as multi:
                multi.write("#!/bin/bash" + "\n"+ "\n")
                multi.write(f"cd therm{st}/" + "\n")
                multi.write("sbatch sub_dftb.sh dftb_in &" + "\n"+ "\n")
                for j in range(st+1, end+1,1):
                        multi.write(f"cd ../therm{j}/"+ "\n")
                        multi.write("sbatch sub_dftb.sh dftb_in &" + "\n" + "\n")
                        multi.write("wait" + "\n")
                        multi.write("echo 'All scripts executed (possibly in parallel).'" + "\n")

        with open(os.path.join("sub_therm_multi.sh"), 'w') as therm_multi:
                therm_multi.write("#!/bin/bash" + "\n"+ "\n")
                therm_multi.write(f"cd therm{st}/" + "\n")
                therm_multi.write("sbatch py_therm.sh &" + "\n"+ "\n")
                for j in range(st+1, end+1,1):
                        therm_multi.write(f"cd ../therm{j}/"+ "\n")
                        therm_multi.write("sbatch py_therm.sh &" + "\n" + "\n")
                        therm_multi.write("wait" + "\n")
                        therm_multi.write("echo 'All scripts executed (possibly in parallel).'" + "\n")