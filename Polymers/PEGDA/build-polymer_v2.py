#!/usr/bin/env python
# coding: utf-8

# In[1]:


# In this code we will optimize a monomer Psi4
# Then we wil calculate RESP charge (https://github.com/cdsgroup/resp)
# Then build the polymer using Pysimm
# we have the initial structure of the monomer from avogadro
# Let's first build the monomer and geometry optimize it 
### Written by Supriyo Naskar, ICGM, Univ. of Montepllier, CNRS, ENSCM ####


# In[2]:


from __future__ import division
import math 
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt 
from pysimm import system, lmps, forcefield, amber
from pysimm.apps.random_walk import random_walk
from pysimm.models.monomers.gaff2.pe import monomer
from pysimm.apps.random_walk import random_walk
from pysimm.apps.random_walk import copolymer
import sys
import psi4
import resp
import os
import pandas as pd


# In[176]:


print ("""# In this code we will optimize a monomer structure using scf/cc-pVDZ \n
# Then we wil calculate ESP and RESP charges\n 
# Then we build the polymer \n
# we have the initial structure of the monomer from avogadro \n 
# You need to have the following modules: numpy, math, mdtraj, pysimm, psi4, resp, pandas \n
# you need to have the follwing files: input_parameters.dat file and {monomername}_initial.xyz \n 
# input_parameters.dat file should have 7 columns \n
#1. Index of head atom of the monomer (index starts from 1) \n
#2. Index of the tail atom of the monomer \n
#3. Hydrogen attached to the head atom of the monomer \n
#4. Hydrogen attached to the tail atom of the monomer \n
#5. Name of the molecule (You should remember that your initial xyz file of the monomer should have a name like {moleculename}_monomer.mol2), also make sure xyz file have two remarks like in begining. \n 
#6. Flag to keep or remove all the unnecessary files dusing the code execution. if 0 it keeps everything, elseif 1 it removes everything \n
#7. Flag to optimize or not, if 1 it will optimize, if 0 it will build the polymer...else not 
# Sometimes, you may not have some parameters in your force-field database, then you have to modify gaff.json to include those. \n
### Written by Supriyo Naskar, ICGM, Univ. of Montepllier, CNRS, ENSCM ### """)


# In[3]:


psi4.set_memory('100 GB')
psi4.set_num_threads(20)


# In[188]:


###define the index (starts from 1) of head/tail and attached H atoms 
i=0
for line in open("input_parameters.dat"):
    li=line.strip()
    if not li.startswith("#"):
        if len(line.rstrip()):
            if  i == 0: 
                head1  = int(line.rstrip())
                i = i+1
            elif  i == 1:
                tail1  = int(line.rstrip())
                i = i+1
            elif  i == 2:
                head1_H  = int(line.rstrip())
                i = i+1
            elif  i == 3:
                tail1_H  = int(line.rstrip())
                i = i+1
            elif  i == 4 : 
                molecule_name = str(line.rstrip())
                i =i+1
            elif  i == 5 :
                remove1 = int(line.rstrip())
                i=i+1
            elif  i == 6 :
                opt = int(line.rstrip())
                i=i+1



# In[5]:

if opt ==1: 

    with open(F'{molecule_name}_initial.xyz') as f: 
        next(f)
        next(f)
        xyz = f.read() 
    PBI = psi4.geometry(xyz)
    
    
    # In[6]:
    
    
    # Set the name of the output file for the initial energy calculation
    # Calculate the initial energy of the molecule using the Hartree-Fock method
    # and the cc-pVDZ basis set and print the output to a file
    psi4.set_output_file(F'{molecule_name}_energy_initial.dat', False)
    psi4.energy('scf/cc-pVDZ',molecule=PBI)
    
    # Set the name of the output file to write the geometry information
    # Print atomic coordinates and interatomic distances to this file
    psi4.set_output_file(F'{molecule_name}_geometry_initial.dat', False)
    PBI.print_out_in_angstrom()
    
    
    # In[7]:
    
    
    # optimize the molecular geometry
    psi4.set_output_file(F'{molecule_name}_geometry_optimization.dat', False)
    psi4.set_options({
        'geom_maxiter': 500,
        'maxiter': 500
    })
    psi4.optimize('scf/cc-pVDZ', molecule=PBI)
    
    # print the optimized atomic coordinates and interatomic distances
    psi4.set_output_file(F'{molecule_name}_geometry_final.dat', False)
    PBI.print_out_in_angstrom()
    PBI.save_xyz_file(F'{molecule_name}_optimized.xyz',1)
    
    
    
    
    with open(F'{molecule_name}_optimized.xyz') as f: 
        next(f)
        next(f)
        xyz2 = f.read() 
    
    
    # In[9]:
    
    
    mol1 = psi4.geometry(xyz2)
    mol1.update_geometry()
    mol1.set_name('optimized')
    
    
    # In[10]:
    
    
    options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
               'VDW_POINT_DENSITY'  : 1.0,
               'RESP_A'             : 0.0005,
               'RESP_B'             : 0.1,
               }
    charges = resp.resp([mol1], options)
    
    
    # In[11]:
    
    
    options['RESP_A'] = 0.001
    options = {'constraint_charge': [[0.0, [head1_H]], [0.0, [tail1_H]]]}
    charges2 = resp.resp([mol1], options)
    
    
    # In[ ]:
    
    
    
    
    
    # In[ ]:
    
    
    
    
    
    # In[12]:
    
    
    ### write the mol2 for Pysimm ###
    traj = pd.read_table(F'{molecule_name}_optimized.xyz', skiprows=2, delim_whitespace=True, names=['atom','x', 'y', 'z'])
    no_atoms = len (traj.x)
    f = open(F'{molecule_name}_monomer.mol2','w+')
    
    bond1 = []
    bond2 = []
    bond3 = []
    bond4 = []
    k = 0 
    for i in range (no_atoms):
        at1 = traj.atom[i]
        for j in range (i+1,no_atoms):
            at2 = traj.atom[j]
            x1=traj.x[i]
            y1=traj.y[i]
            z1=traj.z[i]
            x2=traj.x[j]
            y2=traj.y[j]
            z2=traj.z[j]
            dis= np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
            if dis < 1.9 and (at1 != 'H' and at2 != 'H')   :
                k = k+1
                bond1.append(k) 
                bond2.append(i+1)
                bond3.append(j+1)
                bond4.append(1)
            elif dis < 1.6:
                k = k+1
                bond1.append(k)
                bond2.append(i+1)
                bond3.append(j+1)
                bond4.append(1)
    f.write("@<TRIPOS>MOLECULE\n")
    f.write("generated by SUPRIYO\n")
    f.write("  %d   %d    %d     0     0\n"%(no_atoms,len(bond1),0))
    f.write("SMALL\n")
    f.write("USER_CHARGES\n")
    f.write("****\n")
    f.write("Energy = 0\n")
    f.write("\n@<TRIPOS>ATOM\n")
    resname = 'XXX'

    for i in range (no_atoms):
        at = traj.atom[i]
        x=traj.x[i]
        y=traj.y[i]
        z=traj.z[i]
        f.write("%7d %-4s%14.4f%10.4f%10.4f %-2s%10d %-3s%14.6f\n"  %(i+1,at,x,y,z,at,i+1,resname,charges2[1][i]))
    f.write("@<TRIPOS>BOND\n")
    for i in range (len(bond1)):
        f.write("%5d%6d%6d%3d\n"%(bond1[i],bond2[i],bond3[i],bond4[i]))
    f.write("\n@<TRIPOS>SUBSTRUCTURE\n")
    for i in range (1):
        f.write("1 ****        1 TEMP                        0 ****  **** 0 ROOT \n")   
    f.close()





### now do the Pysimmm Calculation ####
if opt ==0:
    def monomer1(is_capped=True):
        s = system.read_mol2(F'{molecule_name}_monomer.mol2')
        m = s.molecules[1]
        amber.get_forcefield_types_cop(s, types='gaff',f=forcefield.Gaff(),name="PEGDA")
        for b in s.bonds:
            if b.a.bonds.count == 4 and b.b.bonds.count == 4:
                b.order = 5
        s.apply_forcefield(forcefield.Gaff())
        c1 = s.particles[head1]
        c2 = s.particles[tail1]
        c1.linker = 'head'
        c2.linker = 'tail'
        lmps.quick_min(s, min_style='fire', save_input=False )
        s.add_particle_bonding()
        return s
    
    
    pe = monomer1(is_capped=False)
    f = forcefield.Gaff()

    def monomer2(is_capped=True):
        s = system.read_mol2('Head.mol2')
        m = s.molecules[1]
        amber.get_forcefield_types_cop(s, types='gaff',f=forcefield.Gaff(),name="Head")
        for b in s.bonds:
            if b.a.bonds.count == 4 and b.b.bonds.count == 4:
                b.order = 5
        s.apply_forcefield(forcefield.Gaff())
        c1 = s.particles[1]
        c2 = s.particles[2]
        c1.linker = 'head'
        c2.linker = 'tail'
        lmps.quick_min(s, np=16, min_style='fire', save_input=False)
        s.add_particle_bonding()
        return s

    he = monomer2(is_capped=False)
    f = forcefield.Gaff()

    def monomer3(is_capped=True):
        s = system.read_mol2('Tail.mol2')
        m = s.molecules[1]
        amber.get_forcefield_types_cop(s, types='gaff',f=forcefield.Gaff(),name="Tail")
        for b in s.bonds:
            if b.a.bonds.count == 4 and b.b.bonds.count == 4:
                b.order = 5
        s.apply_forcefield(forcefield.Gaff())
        c1 = s.particles[1]
        c2 = s.particles[2]
        c1.linker = 'head'
        c2.linker = 'tail'
        lmps.quick_min(s, np=16, min_style='fire', save_input=False)
        s.add_particle_bonding()
        return s

    te = monomer3(is_capped=False)
    f = forcefield.Gaff()

    
    
    name_polymer = []
    freq_polymer=[]
    for i in range (10,51): 
        k = i + 2
        list1 = []
        list1.append("te")
        for j in range (i):
            list1.append("pe")
        list1.append("he")
        list2 = ', '.join(list1)
        exec(f'polymer_{i} = copolymer([{list2}], {k}, forcefield=forcefield.Gaff(),density = 0.001)')
        exec(f'polymer_{i}.center(what="particles", at=[0.0, 0.0, 0.0], move_both=False)')
        exec(f'sim1 = lmps.Simulation(polymer_{i}, log="steps.log")')
        sim1.add_min(min_style='cg', name='minimization',  etol=1.0e-5, ftol=1.0e-5)
        sim1.add_md(ensemble='nvt', timestep=1, temperature=300, length =5000, new_v=True)
        sim1.run(np=16)
        exec(f'polymer_{i}.write_lammps("polymer_{i}.lmps")')
        exec(f'polymer_{i}.write_pdb("polymer_{i}.pdb")')

if remove1 == 1 : 
    try:
        os.remove ('1_default_grid.dat')
    except OSError:
        pass
    try:
        os.remove ('1_optimized_grid.dat')
    except OSError:
        pass
    try:
        os.remove ('1_default_grid_esp.dat')
    except OSError:
        pass
    try:
        os.remove ('1_optimized_grid_esp.dat')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_AC.AC')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_AC.AC0')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_AM1BCC.AC')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_AM1BCC_PRE.AC')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_BOND_TYPE.AC')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_BOND_TYPE.AC0')
    except OSError:
        pass
    try:
        os.remove ('ANTECHAMBER_AM1BCC_PRE.AC')
    except OSError:
        pass
    try:
        os.remove ('ATOMTYPE.INF')
    except OSError:
        pass
    try:
        os.remove ('log.lammps')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_energy_initial.dat')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_energy_initial.log')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_final.dat')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_final.log')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_initial.dat')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_initial.log')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_optimization.dat')
    except OSError:
        pass
    try:
        os.remove (F'{molecule_name}_geometry_optimization.log')
    except OSError:
        pass
    try:
        os.remove ('pysimm.sim.in')
    except OSError:
        pass
    try:
        os.remove ('pysimm.tmp.ac')
    except OSError:
        pass
    try:
        os.remove ('pysimm.tmp.pdb')
    except OSError:
        pass
    try:
        os.remove ('random_walk.xyz')
    except OSError:
        pass
    try:
        os.remove ('relax.log')
    except OSError:
        pass
    try:
        os.remove ('results.out')
    except OSError:
        pass
    try:
        os.remove ('sqm.in')
    except OSError:
        pass
    try:
        os.remove ('sqm.out')
    except OSError:
        pass
    try:
        os.remove ('sqm.pdb')
    except OSError:
        pass
    try:
        os.remove ('capping_opt.log')
    except OSError:
        pass
    try:
        os.remove ('timer.dat')
    except OSError:
        pass






