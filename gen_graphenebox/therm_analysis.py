import os
import numpy as np
import math

def extract_md():
    md_step = []
    pot_energy = []
    kin_energy =[]
    tot_energy = []
    md_temp =[]
    md_path = os.path.join("md.out")
    with open(md_path, 'r') as md_out:
        data= md_out.readlines()
        for indx,line in enumerate(data):
            text = line.split(':')
            if 'MD step' in text:
                md_step.append(int(text[1]))
            if 'Potential Energy' in text:
                pot_energy.append(float(text[1].split()[0]))
            if 'MD Kinetic Energy' in text:
                kin_energy.append(float(text[1].split()[0]))
            if 'Total MD Energy' in text:
                tot_energy.append(float(text[1].split()[0]))
            if 'MD Temperature' in text:
                md_temp.append(float(text[1].split()[2]))
    return md_step, pot_energy, kin_energy, tot_energy, md_temp


def md_analysis(n, start, end):
    nmovingatoms = 25 + 25 +(3*n)+4 #for relaxed dual cation system
#    nmovingatoms = (3*n)+2 #for single cation system
#    nmovingatoms = (3*n)+4 #for double cation system
#    nmovingatoms = 19 + (3*n)+2 #for relaxed single cation system
#    nmovingatoms = 50 + (3*n)+4 #for relaxed double cation system
    step, pot, kin, tot, temp =extract_md()
    st=int(start/10)
    en= int((end+10)/10)

    pmean = np.mean(pot[st:en])
    pstd = np.std(pot[st:en])

    kmean = np.mean(kin[st:en])
    kstd = np.std(kin[st:en])

    tmean = np.mean(tot[st:en])
    tstd = np.std(tot[st:en])

    with open(os.path.join("therm.dat"), 'w') as therm_out:
        therm_out.write('temperature from mean kinetic energy          = ' + str(kmean*2*315777/(nmovingatoms*3)) +'K' + "\n")
        therm_out.write('thermal fluctuations of the total energy       = '+ str(1*100/np.sqrt(3*nmovingatoms)) + '%' + "\n")
        therm_out.write('total energy fluctuations compared to 3N*KT   = ' + str(tstd*100/(kmean*2)) + '%' + "\n")
        therm_out.write('kinetic energy fluctuations compared to 3N*KT = ' + str(kstd*100/(kmean*2)) + '%')