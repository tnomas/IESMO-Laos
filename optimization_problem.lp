\* Source Pyomo model name=unknown *\

min 
cost:
+400000 Ppv
+1000000 Pwind
+0.001 UsageWater(1)
+0.001 UsageWater(2)
+0.001 UsageWater(3)
+0.001 UsageWater(4)
+0.001 UsageWater(5)

s.t.

c_l_Energy(1)_:
+0.20000000000000001 Ppv
+0.5 Pwind
+0.001 UsageWater(1)
>= 500

c_l_Energy(2)_:
+0.29999999999999999 Ppv
+0.29999999999999999 Pwind
+0.001 UsageWater(2)
>= 400

c_l_Energy(3)_:
+0.80000000000000004 Ppv
+0.10000000000000001 Pwind
+0.001 UsageWater(3)
>= 300

c_l_Energy(4)_:
+0.40000000000000002 Ppv
+0.59999999999999998 Pwind
+0.001 UsageWater(4)
>= 450

c_l_Energy(5)_:
+0.02 Ppv
+0.20000000000000001 Pwind
+0.001 UsageWater(5)
>= 310

c_u_MaxWaterPower(1)_:
+0.001 UsageWater(1)
<= 480

c_u_MaxWaterPower(2)_:
+0.001 UsageWater(2)
<= 480

c_u_MaxWaterPower(3)_:
+0.001 UsageWater(3)
<= 480

c_u_MaxWaterPower(4)_:
+0.001 UsageWater(4)
<= 480

c_u_MaxWaterPower(5)_:
+0.001 UsageWater(5)
<= 480

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
   0 <= UsageWater(1) <= +inf
   0 <= UsageWater(2) <= +inf
   0 <= UsageWater(3) <= +inf
   0 <= UsageWater(4) <= +inf
   0 <= UsageWater(5) <= +inf
end
