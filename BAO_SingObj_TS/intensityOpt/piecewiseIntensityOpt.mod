option solver ipopt; 
options ipopt_options "linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001"; 
			#option solver gurobi_ampl; 
param R1; 		#number of voxels of PTVHD
param bmlt; 		#total number of beamlets
param bmlt_opt; 	#total number of beamlets in the optimal solution
var x {1 .. bmlt} >= 0, <=400, default 1;
param ddmPTVHD{1 .. R1, 1 .. bmlt};
param ddmPTVHD_opt{1 .. R1, 1 .. bmlt_opt};
param x_opt{1 .. bmlt_opt};
var pVar {i in 1..R1} >=0;
var nVar {i in 1..R1} >=0;

minimize intensityDiff: sum {i in 1..R1} (pVar[i] + nVar[i]);

s.t. 

beamletsEq_1 {j in 71..bmlt}: x_opt[j] = x[j];
absConstr1 {i in 1 .. R1}: (sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])-(sum {j in 1..bmlt_opt} x_opt[j]*ddmPTVHD_opt[i,j]) = pVar[i] - nVar[i];
#absConstr2 {i in 1 .. R1}: y_2[i] = (-(sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])+(sum {j in 1..bmlt_opt} x_opt[j]*ddmPTVHD_opt[i,j]) )* (1-z[i]);

