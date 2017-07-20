\* Source Pyomo model name=unknown *\

min 
cost:
+50000 Ppv
+50000 Pwind
<<<<<<< HEAD
+9640000 ONE_VAR_CONSTANT
=======
+9000000 ONE_VAR_CONSTANT
>>>>>>> 9bc24b276eabb77a1a7ad11b040c5777afe56346

s.t.

c_l_EnergyDemand(1)_:
+4.0048400000000001e-05 Pwind
<<<<<<< HEAD
>= -168.39911746707841

c_l_EnergyDemand(2)_:
+7.6765600000000002e-05 Pwind
>= -168.21548409103679

c_l_EnergyDemand(3)_:
+4.4193900000000004e-05 Pwind
>= -167.7563934771936

c_l_EnergyDemand(4)_:
+1.6619000000000001e-05 Pwind
>= -167.02185279928801

c_l_EnergyDemand(5)_:
+1.3990999999999999e-05 Pwind
>= -166.195473912144
=======
>= 0.60088253292160054

c_l_EnergyDemand(2)_:
+7.6765600000000002e-05 Pwind
>= 0.78451590896320056

c_l_EnergyDemand(3)_:
+4.4193900000000004e-05 Pwind
>= 1.2436065228064006

c_l_EnergyDemand(4)_:
+1.6619000000000001e-05 Pwind
>= 1.9781472007120016

c_l_EnergyDemand(5)_:
+1.3990999999999999e-05 Pwind
>= 2.8045260878559999
>>>>>>> 9bc24b276eabb77a1a7ad11b040c5777afe56346

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
end
