\* Source Pyomo model name=unknown *\

min 
cost:
+26666.666666666668 Ppv
+50000 Pwind

s.t.

c_l_demand(1)_:
+0.20000000000000001 Ppv
+0.5 Pwind
>= 5

c_l_demand(2)_:
+0.29999999999999999 Ppv
+0.29999999999999999 Pwind
>= 9000

c_l_demand(3)_:
+0.80000000000000004 Ppv
+0.10000000000000001 Pwind
>= 20

c_l_demand(4)_:
+0.40000000000000002 Ppv
+0.59999999999999998 Pwind
>= 21

c_l_demand(5)_:
+0.20000000000000001 Pwind
>= 300

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
end
