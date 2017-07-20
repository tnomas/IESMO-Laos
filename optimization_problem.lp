\* Source Pyomo model name=unknown *\

min 
cost:
+0.01 PowerGeneratingWater(1)
+0.01 PowerGeneratingWater(10)
+0.01 PowerGeneratingWater(11)
+0.01 PowerGeneratingWater(12)
+0.01 PowerGeneratingWater(13)
+0.01 PowerGeneratingWater(14)
+0.01 PowerGeneratingWater(15)
+0.01 PowerGeneratingWater(16)
+0.01 PowerGeneratingWater(17)
+0.01 PowerGeneratingWater(18)
+0.01 PowerGeneratingWater(19)
+0.01 PowerGeneratingWater(2)
+0.01 PowerGeneratingWater(20)
+0.01 PowerGeneratingWater(21)
+0.01 PowerGeneratingWater(22)
+0.01 PowerGeneratingWater(23)
+0.01 PowerGeneratingWater(24)
+0.01 PowerGeneratingWater(3)
+0.01 PowerGeneratingWater(4)
+0.01 PowerGeneratingWater(5)
+0.01 PowerGeneratingWater(6)
+0.01 PowerGeneratingWater(7)
+0.01 PowerGeneratingWater(8)
+0.01 PowerGeneratingWater(9)
+3333.3333333333335 Ppv
+50000 Pwind

s.t.

c_l_EnergyDemand(1)_:
+0.0001 PowerGeneratingWater(1)
+4.0048400000000001e-05 Pwind
>= 5.6008825329216005

c_l_EnergyDemand(2)_:
+0.0001 PowerGeneratingWater(2)
+7.6765600000000002e-05 Pwind
>= 5.7845159089632006

c_l_EnergyDemand(3)_:
+0.0001 PowerGeneratingWater(3)
+4.4193900000000004e-05 Pwind
>= 6.2436065228064006

c_l_EnergyDemand(4)_:
+0.0001 PowerGeneratingWater(4)
+1.6619000000000001e-05 Pwind
>= 6.9781472007120016

c_l_EnergyDemand(5)_:
+0.0001 PowerGeneratingWater(5)
+1.3990999999999999e-05 Pwind
>= 7.8045260878559999

c_l_EnergyDemand(6)_:
+0.0001 PowerGeneratingWater(6)
+1.6619000000000001e-05 Pwind
>= 8.3553975210240008

c_l_EnergyDemand(7)_:
+0.0001 PowerGeneratingWater(7)
+0.00013151852264499999 Ppv
+1.92469e-05 Pwind
>= 8.6308691063040008

c_l_EnergyDemand(8)_:
+0.0001 PowerGeneratingWater(8)
+0.00026086318541099999 Ppv
+3.2571699999999999e-05 Pwind
>= 8.9981645533440009

c_l_EnergyDemand(9)_:
+0.0001 PowerGeneratingWater(9)
+0.00037173003921099997 Ppv
+4.8339399999999993e-05 Pwind
>= 8.7226929680640009

c_l_EnergyDemand(10)_:
+0.0001 PowerGeneratingWater(10)
+0.00045270631968199993 Ppv
+8.7906599999999996e-05 Pwind
>= 8.6308691063040008

c_l_EnergyDemand(11)_:
+0.0001 PowerGeneratingWater(11)
+0.00049129233252400003 Ppv
+0.0001569365 Pwind
>= 8.6308691063040008

c_l_EnergyDemand(12)_:
+0.0001 PowerGeneratingWater(12)
+0.00048857500767600002 Ppv
+0.0002622023 Pwind
>= 9.2735644012320009

c_l_EnergyDemand(13)_:
+0.0001 PowerGeneratingWater(13)
+0.00044401088016900002 Ppv
+0.00033071400000000002 Pwind
>= 8.9981645533440009

c_l_EnergyDemand(14)_:
+0.0001 PowerGeneratingWater(14)
+0.00036140420478800005 Ppv
+0.00031886970000000002 Pwind
>= 8.3553975210240008

c_l_EnergyDemand(15)_:
+0.0001 PowerGeneratingWater(15)
+0.00024836349111000003 Ppv
+0.00021067979999999999 Pwind
>= 8.2635736592640008

c_l_EnergyDemand(16)_:
+0.0001 PowerGeneratingWater(16)
+0.000116301503496 Ppv
+9.4013799999999998e-05 Pwind
>= 8.2635736592640008

c_l_EnergyDemand(17)_:
+0.0001 PowerGeneratingWater(17)
+1.03258344225277e-05 Ppv
+4.4193900000000004e-05 Pwind
>= 9.0899166777120008

c_l_EnergyDemand(18)_:
+0.0001 PowerGeneratingWater(18)
+4.0048400000000001e-05 Pwind
>= 9.916331433552001

c_l_EnergyDemand(19)_:
+0.0001 PowerGeneratingWater(19)
+3.2571699999999999e-05 Pwind
>= 10.926322175520001

c_l_EnergyDemand(20)_:
+0.0001 PowerGeneratingWater(20)
+4.0048400000000001e-05 Pwind
>= 10.559026728480001

c_l_EnergyDemand(21)_:
+0.0001 PowerGeneratingWater(21)
+5.6630399999999999e-05 Pwind
>= 9.7326837100320009

c_l_EnergyDemand(22)_:
+0.0001 PowerGeneratingWater(22)
+6.6698000000000004e-05 Pwind
>= 8.7226929680640009

c_l_EnergyDemand(23)_:
+0.0001 PowerGeneratingWater(23)
+6.1664199999999998e-05 Pwind
>= 7.8963499496159999

c_l_EnergyDemand(24)_:
+0.0001 PowerGeneratingWater(24)
+5.2484899999999996e-05 Pwind
>= 5.8763325969840006

c_u_MaxWaterPower(1)_:
+0.0001 PowerGeneratingWater(1)
<= 5

c_u_MaxWaterPower(2)_:
+0.0001 PowerGeneratingWater(2)
<= 5

c_u_MaxWaterPower(3)_:
+0.0001 PowerGeneratingWater(3)
<= 5

c_u_MaxWaterPower(4)_:
+0.0001 PowerGeneratingWater(4)
<= 5

c_u_MaxWaterPower(5)_:
+0.0001 PowerGeneratingWater(5)
<= 5

c_u_MaxWaterPower(6)_:
+0.0001 PowerGeneratingWater(6)
<= 5

c_u_MaxWaterPower(7)_:
+0.0001 PowerGeneratingWater(7)
<= 5

c_u_MaxWaterPower(8)_:
+0.0001 PowerGeneratingWater(8)
<= 5

c_u_MaxWaterPower(9)_:
+0.0001 PowerGeneratingWater(9)
<= 5

c_u_MaxWaterPower(10)_:
+0.0001 PowerGeneratingWater(10)
<= 5

c_u_MaxWaterPower(11)_:
+0.0001 PowerGeneratingWater(11)
<= 5

c_u_MaxWaterPower(12)_:
+0.0001 PowerGeneratingWater(12)
<= 5

c_u_MaxWaterPower(13)_:
+0.0001 PowerGeneratingWater(13)
<= 5

c_u_MaxWaterPower(14)_:
+0.0001 PowerGeneratingWater(14)
<= 5

c_u_MaxWaterPower(15)_:
+0.0001 PowerGeneratingWater(15)
<= 5

c_u_MaxWaterPower(16)_:
+0.0001 PowerGeneratingWater(16)
<= 5

c_u_MaxWaterPower(17)_:
+0.0001 PowerGeneratingWater(17)
<= 5

c_u_MaxWaterPower(18)_:
+0.0001 PowerGeneratingWater(18)
<= 5

c_u_MaxWaterPower(19)_:
+0.0001 PowerGeneratingWater(19)
<= 5

c_u_MaxWaterPower(20)_:
+0.0001 PowerGeneratingWater(20)
<= 5

c_u_MaxWaterPower(21)_:
+0.0001 PowerGeneratingWater(21)
<= 5

c_u_MaxWaterPower(22)_:
+0.0001 PowerGeneratingWater(22)
<= 5

c_u_MaxWaterPower(23)_:
+0.0001 PowerGeneratingWater(23)
<= 5

c_u_MaxWaterPower(24)_:
+0.0001 PowerGeneratingWater(24)
<= 5

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds
   0 <= Pwind <= +inf
   0 <= Ppv <= +inf
   0 <= PowerGeneratingWater(1) <= +inf
   0 <= PowerGeneratingWater(2) <= +inf
   0 <= PowerGeneratingWater(3) <= +inf
   0 <= PowerGeneratingWater(4) <= +inf
   0 <= PowerGeneratingWater(5) <= +inf
   0 <= PowerGeneratingWater(6) <= +inf
   0 <= PowerGeneratingWater(7) <= +inf
   0 <= PowerGeneratingWater(8) <= +inf
   0 <= PowerGeneratingWater(9) <= +inf
   0 <= PowerGeneratingWater(10) <= +inf
   0 <= PowerGeneratingWater(11) <= +inf
   0 <= PowerGeneratingWater(12) <= +inf
   0 <= PowerGeneratingWater(13) <= +inf
   0 <= PowerGeneratingWater(14) <= +inf
   0 <= PowerGeneratingWater(15) <= +inf
   0 <= PowerGeneratingWater(16) <= +inf
   0 <= PowerGeneratingWater(17) <= +inf
   0 <= PowerGeneratingWater(18) <= +inf
   0 <= PowerGeneratingWater(19) <= +inf
   0 <= PowerGeneratingWater(20) <= +inf
   0 <= PowerGeneratingWater(21) <= +inf
   0 <= PowerGeneratingWater(22) <= +inf
   0 <= PowerGeneratingWater(23) <= +inf
   0 <= PowerGeneratingWater(24) <= +inf
end
