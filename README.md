# Targetted-Recombination
Code for 1) building the constraints matrix using an Rccp function and 2) solving the integer program with the constraints to determine the optimal recombination points in an F1 to maximize gain.
 
The inputs are 	1) Matrix of phased F1 marker genotypes
			i.e. 1	0
			     0	1
			     1  1
 		2) Vector of marker-effect estimates
			i.e. 0.01
			    -0.03
			     0.008
		3) A genetic map associating markers with positions on linkage groups.
			i.e.
			    Marker	Chr
			     1		1
			     2		1
			     3		2
		
