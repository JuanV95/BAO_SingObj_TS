option solver ipopt; 
options ipopt_options "linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001"; 
#option solver gurobi_ampl; 
param R1; #number of voxels of PTVHD
#param R2; #number of voxels of RECTUM
#param R3; #number of voxels of BLADDER
param bmlt; 	#total number of beamlets
param bmlt_opt; 	#total number of beamlets in the optimal solution
var x {1 .. bmlt} >= 0, <=400, default 1;
#var y_1 {1 .. R1} >= 0, <=400;
#var y_2 {1 .. R1} >= 0, <=400;  
#var z {1 .. R1} binary;  
param ddmPTVHD{1 .. R1, 1 .. bmlt};
#param ddmRECTUM{1 .. R2, 1 .. bmlt};
#param ddmBLADDER{1 .. R3, 1 .. bmlt};
param ddmPTVHD_opt{1 .. R1, 1 .. bmlt_opt};
#param ddmRECTUM_opt{1 .. R2, 1 .. bmlt_opt};
#param ddmBLADDER_opt{1 .. R3, 1 .. bmlt_opt};
param x_opt{1 .. bmlt_opt};

minimize intensityDiff: sum {i in 1..R1} abs((sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])-(sum {j in 1..bmlt_opt} x_opt[j]*ddmPTVHD_opt[i,j]));

s.t. 

beamletsEq_1 {j in 71..bmlt}: x_opt[j] = x[j];

#absConstr1 {i in 1 .. R1}: (sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])-(sum {j in 1..bmlt_opt} x_opt[j]*ddmPTVHD_opt[i,j]) = pVar[i] + nVar[i];
#absConstr2 {i in 1 .. R1}: y_2[i] = (-(sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])+(sum {j in 1..bmlt_opt} x_opt[j]*ddmPTVHD_opt[i,j]) )* (1-z[i]);

