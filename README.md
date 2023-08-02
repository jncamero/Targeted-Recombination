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
			
R environment details
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8      
 [2] LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8      
 [8] LC_NAME=C                 
 [9] LC_ADDRESS=C              
[10] LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices
[4] utils     datasets  methods  
[7] base     

other attached packages:
[1] Rcpp_1.0.11 Rglpk_0.6-5 slam_0.1-48

loaded via a namespace (and not attached):
[1] compiler_4.1.3 tools_4.1.3 

*GLPK utilities to install before installing Rglpk
apt-get install glpk-utils libglpk-dev glpk-doc
