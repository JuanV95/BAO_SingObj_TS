option solver ipopt; 
options ipopt_options "linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001"; 
param R1; #number of voxels of PTVHD
param R2; #number of voxels of RECTUM
param R3; #number of voxels of BLADDER
param bmlt; 	#number of beamlets 
param ddmPTVHD{1 .. R1, 1 .. bmlt};
param ddmRECTUM{1 .. R2, 1 .. bmlt};
param ddmBLADDER{1 .. R3, 1 .. bmlt};
param UB1;
param UB2;
param UB3;
param t;
param epsilon;
param OAR_targetUB;
var x {1 .. bmlt} >= 0, <=400, default 1; 
param a{1 .. 4}; 
param v{1 .. 4}; 
param EUD0{1 .. 4}; 
minimize Total_Cost: - log((1+(((1/R2)*(sum {i in 1..R2} (sum {j in 1..bmlt} x[j]*ddmRECTUM[i,j])^a[2]))^(1/a[2])/EUD0[2])^v[2])^-1)- log((1+(((1/R3)*(sum {i in 1..R3} (sum {j in 1..bmlt} x[j]*ddmBLADDER[i,j])^a[3]))^(1/a[3])/EUD0[3])^v[3])^-1);
 s.t. 
equalityTarget: 	((1/R1)*(sum {i in 1..R1} (sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])^a[1]))^(1/a[1]) = t; 
constraintOAR_Target: 	((1/R1)*(sum {i in 1..R1} (sum {j in 1..bmlt} x[j]*ddmPTVHD[i,j])^a[4]))^(1/a[4]) <=OAR_targetUB;
#OAR_UB2: 	((1/R2)*(sum {i in 1..R2} (sum {j in 1..bmlt} x[j]*ddmRECTUM[i,j])^a[2]))^(1/a[2]) <= UB2; 
#OAR_UB3: 	((1/R3)*(sum {i in 1..R3} (sum {j in 1..bmlt} x[j]*ddmBLADDER[i,j])^a[3]))^(1/a[3]) <= UB3; 
