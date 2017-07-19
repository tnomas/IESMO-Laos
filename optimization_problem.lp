\* Source Pyomo model name=unknown *\

min 
cost:
+50000 Ppv
+50000 Pwind

s.t.

c_l_EnergyDemand(1)_:
+0.5 Ppv
+4.0048400000000001e-05 Pwind
>= 7.80748e-05

c_l_EnergyDemand(2)_:
+0.29999999999999999 Ppv
+7.6765600000000002e-05 Pwind
>= 8.0634600000000002e-05

c_l_EnergyDemand(3)_:
+0.80000000000000004 Ppv
+4.4193900000000004e-05 Pwind
>= 8.7034199999999994e-05

c_l_EnergyDemand(4)_:
+0.40000000000000002 Ppv
+1.6619000000000001e-05 Pwind
>= 9.727350000000001e-05

c_l_EnergyDemand(5)_:
+1.3990999999999999e-05 Pwind
>= 0.00010879299999999999

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
end
