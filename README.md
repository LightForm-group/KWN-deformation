
@author Madeleine Bignon, University of Manchester  
@author Pratheek Shanthraj, University of Manchester  
@author Samuel Engel, University of Manchester  

KWN precipitation model including the effect of deformation via excess vacancy concentration  
The code runs the model described in ref [3]


References:   
 [1] Robson, J. D. (2020). Metallurgical and Materials Transactions A: Physical Metallurgy and Materials Science, 51(10), 5401–5413. https://doi.org/10.1007/s11661-020-05960-5  
 [2] Deschamps, A., Livet, F., & Bréchet, Y. (1998). . Acta Materialia, 47(1), 281–292. https://doi.org/10.1016/S1359-6454(98)00293-6  
 [3] Bignon, M., Shanthraj, P., & Robson, J. D. (2022). Acta Materialia, 234, 118036. https://doi.org/10.1016/J.ACTAMAT.2022.118036   
 [4] Deschamps, A., & De Geuser, F. (2011). Journal of Applied Crystallography, 44(2), 343–352. https://doi.org/10.1107/S0021889811003049  
 [5] Perez, M. (2005). Scripta Materialia, 52(8), 709–712. https://doi.org/10.1016/j.scriptamat.2004.12.026  
 [6] Nicolas, M., & Deschamps, A. (2003).  Acta Materialia, 51(20), 6077–6094. https://doi.org/10.1016/S1359-6454(03)00429-4  
  
-------------------------------------------------------------------------------------------   
KWN model with external deformation considered via the amount of excess vacancies  
-------------------------------------------------------------------------------------------   
This model is described in details in ref [3]  
The program allows to calculate the evolution of a precipitate population in static conditions or under deformation  
A classical KWN framework is used to calculate the precipitation kinetics and the effect of deformation is considered through enhanced diffusivity due to the production of excess vacancies. 
The effect of deformation on excess vacancies is considered via a phenomenological model.    


The code has been written in MacOS 11.6 and compiled with ifort. 
If used with Windows, the fortran file (KWN_with_deformation.f90) should be recompiled with ifort. 

  
To run the code :  
1. Fill or modify the "input.yaml" file with input corresponding to the model described in ref. [3].
2. Run ./run_kwn.sh in a command line.  
3. The outputs are written in textfiles and can be visualised using the attached Jupyter notebook.  


Some examples of input files canb= be found in the "example" directory: 
- with deformation  
   - for two solute elements
   - for binary alloys 
- without deformation 
   - or two solute elements
   - for binary alloys 
   
