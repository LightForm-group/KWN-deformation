import numpy as np
import pandas as pd
import csv
import math
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import quad
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import os

test_folder='test_nickel_dissolution'

initial_guess_interfacial_energy=0.09

temperature_max=[1090,1110]
strain_rate=[0]

prefix='{f}/results/'.format(f=test_folder)

initial_precipitate_distribution=[]

suffix=[]






for j in range(len(strain_rate)):
    for i  in range(len(temperature_max)):
        #the suffix of the result files depends on the temperature and strain rate considered
        suffix.append("{temperature}C_strain_rate{str:3.3E}.txt".format(temperature=temperature_max[i], str=strain_rate[j]))




#find the line where the interfacial energy is defined
def run_model(interfacial_energy, temperature):
    replaced_content=""
    word = 'gamma_coherent' #name of interfacial energy in namelist.input
    with open(r'{f}/namelist.input'.format(f=test_folder), 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
    
        for line in lines:
            
            # check if string present on a current line
            if line.find(word) != -1:
               # print(word, 'string exists in file')
               # print('Line Number:', lines.index(line))
               # print('Line:', line)
                line = line.strip()
                new_line=line.replace(line, 'gamma_coherent = {num}'.format(num=interfacial_energy))
                
            elif line.find('Temperature') != -1:
                line = line.strip()
                new_line=line.replace(line, 'Temperature = {num}'.format(num=temperature+273))
            
            else:
                line = line.strip()
                new_line=line
            
            replaced_content = replaced_content + new_line + "\n"
    #change the input file to replace the interfacial energy by desired value
    write_file=open('{f}/namelist.input'.format(f=test_folder),"w")
    write_file.write(replaced_content)
    write_file.close()

    #run the KWN model
    os.system('./run_KWN.sh {f}'.format(f=test_folder))
    
    #copy the results in appropriate location 
    path = "{dir}/optimisation/{it}/".format( dir =test_folder, it=iteration)
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    #print(isExist)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("The new directory is created!")
    suffix="{temperature}C_strain_rate{str:3.3E}.txt".format(temperature=temperature, str=strain_rate[j])

    file="{pre}KWN_parameters_{end}".format(pre=prefix,end=suffix)
    os.system('cp {file} {path}/'.format(file=file, path=path))
    file= "{pre}kinetics_data_{suffix}".format(pre=prefix,suffix=suffix)
    os.system('cp {file} {path}/'.format(file=file, path=path))



def calculate_residual(iteration, temperature):
    path = "{dir}/optimisation/{it}/".format( dir =test_folder, it=iteration)
    suffix="{temperature}C_strain_rate{str:3.3E}.txt".format(temperature=temperature, str=strain_rate[j])
    file_name= "{pre}kinetics_data_{suffix}".format(dir=test_folder, pre=path,suffix=suffix)
        
    # print(file_name)
    data=[]
    data=np.genfromtxt(file_name, skip_header=2)
    time=data[:,0]
    mean_radius=data[:,1] #nm       
    precipitate_density=data[:,2] #/micron^3
    vf=data[:,3] #

    name_experimental_file='{dir}/experimental/time_vf_{temp}.txt'.format(dir=test_folder, temp=temperature)
    data_exp=np.genfromtxt(name_experimental_file)
    vf_exp=data_exp[:,1]
    print('vf exp', vf)


    f=interpolate.interp1d(time, vf)
    x_new=data_exp[:,0]

    #remove data out of calculation range to compute the residual
    vf_exp=vf_exp[x_new>min(time)]
    vf_exp=vf_exp[x_new<max(time)]
    
    x_new=x_new[x_new>min(time)]
    x_new=x_new[x_new<max(time)]

    #interpolate the calculated values on the experimental values 
    print(x_new)
    ynew=f(x_new)
    #print(vf_exp-ynew)
    residual=[x_new, (vf_exp/100-ynew)/vf_exp*100, ynew]
   # print(residual)
    return residual

    #plt.plot(data_exp[:,0], data_exp[:,1]/100, 'x',color= color[i], label='T={temperature}°C'.format(temperature=temperature_max[i]))


t=[]
res=[]
calculated_vf=[]

#first run the model with the initial guess
iteration=0

for temperature in temperature_max:

    run_model(initial_guess_interfacial_energy, temperature)

    #calculate the residual
    residual=calculate_residual(iteration, temperature)
    print('residual:', residual)
    #t is the time
    t=np.concatenate((t,residual[0]))
    #res is the difference between calculated and experimental volume fraction
    res=np.concatenate((res,residual[1]))



#now randomly vary the parameter of interest
interfacial_energy=initial_guess_interfacial_energy
delta=0.01*interfacial_energy #arbitrary factor


for iteration in range(5):
    iteration=iteration+1
    temp_res=res
    temp_interfacial_energy=interfacial_energy   
    t=[]
    res=[]
    calculated_vf=[]
    for temperature in temperature_max:
        interfacial_energy=interfacial_energy+delta
        print('Interfacial energy:', interfacial_energy)
        run_model(interfacial_energy, temperature)

        #calculate the residual
        residual=calculate_residual(iteration, temperature)
        print('residual:', residual)
        #t is the time
        t=np.concatenate((t,residual[0]))
        #res is the difference between calculated and experimental volume fraction
        res=np.concatenate((res,residual[1]))
    
    from numpy.linalg import inv
    jacobian=(res-temp_res)/(interfacial_energy-temp_interfacial_energy)
    hessian=(np.matmul(jacobian.T, jacobian))
    delta=-1/hessian*np.matmul(jacobian.T, res)
    