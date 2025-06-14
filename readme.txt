This folder contains all the codes for running the surface current analysis
1. The efield_to_current_top and efield_to_Current_bottom codes are used to convert the sampled nearfield data into currents. Before running the code, ensure the following
	1. The file paths are changed in the code.
	2. The mathematics used in the code takes the current plane and the sampled near field origin to be at 0,0,0. But in some simulations, this might not be the case. So, ensure that the x,y,z are shifted to origin. This can be done by just adding or subtracting the offset. This is also mentioned as  comment in the code.
	3. For calculating the equivalent current values in both top and bottom planes, the z values of the plane should be mentioned correctly. This is also given as a comment in the code. So please ensure to make this change before starting the code.
	

2. The current to farfield codes converts the currents to farfield pattern. Please ensure the paths are corrected before running the code.

3. Jcurrentpca.m code is used for performing the PCA analysis. This code requires the eigenowncode function to run. Ensure that they are in same folder before running. 

4.Plotsandanalysis_jcurrents code contains all codes for performing different analysis on the data such as the totalpower calculation, FFT. Each part of the code is commented for reference in the code itself. So carefully read the codes and use the required snippets on the code.

5.chromaticity1.m is the code to calculate the chromaticity. This function is called inside the plots and analysis code. Make sure to have this in same folder as that file.