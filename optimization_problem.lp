\* Source Pyomo model name=unknown *\

min 
cost:
+26666.666666666668 Ppv
+50000 Pwind
+1000 UsageWater(1)
+1000 UsageWater(2)
+1000 UsageWater(3)
+1000 UsageWater(4)
+1000 UsageWater(5)

s.t.

c_l_Energy(1)_:
+0.20000000000000001 Ppv
+0.5 Pwind
+100 UsageWater(1)
>= 5

c_l_Energy(2)_:
+0.29999999999999999 Ppv
+0.29999999999999999 Pwind
+100 UsageWater(2)
>= 9000

c_l_Energy(3)_:
+0.80000000000000004 Ppv
+0.10000000000000001 Pwind
+100 UsageWater(3)
>= 20

c_l_Energy(4)_:
+0.40000000000000002 Ppv
+0.59999999999999998 Pwind
+100 UsageWater(4)
>= 21

c_l_Energy(5)_:
+0.20000000000000001 Pwind
+100 UsageWater(5)
>= 300

c_u_Drinkwater(1)_:
+0 ONE_VAR_CONSTANT
<= 4500

c_u_Drinkwater(2)_:
+0 ONE_VAR_CONSTANT
<= 4500

c_u_Drinkwater(3)_:
+0 ONE_VAR_CONSTANT
<= 4500

c_u_Drinkwater(4)_:
+0 ONE_VAR_CONSTANT
<= 4500

c_u_Drinkwater(5)_:
+0 ONE_VAR_CONSTANT
<= 4500

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
   0 <= UsageWater(1) <= 5000
   0 <= UsageWater(2) <= 5000
   0 <= UsageWater(3) <= 5000
   0 <= UsageWater(4) <= 5000
   0 <= UsageWater(5) <= 5000
end
