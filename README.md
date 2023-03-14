# Shifted-SPOD

This repository contains the scripts necessary to calculate solutions of the complex Ginzburg-Landau equation and perform analysis using the SPOD and shifted SPOD algorithms.
Reference: doi:10.1017/jfm.2022.790

Requirements:
- MATLAB with Signal Processing Toolbox because of fir2 function. This function can be commented out, if necessary.

Steps to run the case:

1. Execute the run_GL_solution.m script to compute the Ginzburg-Landau solution. Two mat files will be generated.
2. Execute both run_GL_SPOD.m and run_GL_ShiftedSPOD.m to obtain SPOD modes computed with standard and shifted algorithms, respectively. Each script generates one mat file.
3. The script run_GL_Resolvent.m computes and plots resolvent response modes.
4. The script spacetime_xcorr.m plots space-time correlations with respect to the point of interest x_p. 
