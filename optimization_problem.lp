\* Source Pyomo model name=unknown *\

min 
cost:
+37120 Ppv
+77350 Pwind

s.t.

c_l_EnergyDemand(1)_:
+0.0001 PowerGeneratingWater(1)
+0.072133333333332994 Pwind
>= 5.6008825329216005

c_l_EnergyDemand(2)_:
+0.0001 PowerGeneratingWater(2)
+0.13826666666666698 Pwind
>= 5.7845159089632006

c_l_EnergyDemand(3)_:
+0.0001 PowerGeneratingWater(3)
+0.079600000000000004 Pwind
>= 6.2436065228064006

c_l_EnergyDemand(4)_:
+0.0001 PowerGeneratingWater(4)
+0.029933333333332996 Pwind
>= 6.9781472007120016

c_l_EnergyDemand(5)_:
+0.0001 PowerGeneratingWater(5)
+0.0252 Pwind
>= 7.8045260878559999

c_l_EnergyDemand(6)_:
+0.0001 PowerGeneratingWater(6)
+0.029933333333332996 Pwind
>= 8.3553975210240008

c_l_EnergyDemand(7)_:
+0.0001 PowerGeneratingWater(7)
+0.12745656 Ppv
+0.034666666666666998 Pwind
>= 8.6308691063040008

c_l_EnergyDemand(8)_:
+0.0001 PowerGeneratingWater(8)
+0.25280639999999999 Ppv
+0.058666666666666999 Pwind
>= 8.9981645533440009

c_l_EnergyDemand(9)_:
+0.0001 PowerGeneratingWater(9)
+0.36024911999999998 Ppv
+0.087066666666667014 Pwind
>= 8.7226929680640009

c_l_EnergyDemand(10)_:
+0.0001 PowerGeneratingWater(10)
+0.43872443999999999 Ppv
+0.15833333333333299 Pwind
>= 8.6308691063040008

c_u_MaxWaterPower(1)_:
+0.0001 PowerGeneratingWater(1)
<= 100

c_u_MaxWaterPower(2)_:
+0.0001 PowerGeneratingWater(2)
<= 100

c_u_MaxWaterPower(3)_:
+0.0001 PowerGeneratingWater(3)
<= 100

c_u_MaxWaterPower(4)_:
+0.0001 PowerGeneratingWater(4)
<= 100

c_u_MaxWaterPower(5)_:
+0.0001 PowerGeneratingWater(5)
<= 100

c_u_MaxWaterPower(6)_:
+0.0001 PowerGeneratingWater(6)
<= 100

c_u_MaxWaterPower(7)_:
+0.0001 PowerGeneratingWater(7)
<= 100

c_u_MaxWaterPower(8)_:
+0.0001 PowerGeneratingWater(8)
<= 100

c_u_MaxWaterPower(9)_:
+0.0001 PowerGeneratingWater(9)
<= 100

c_u_MaxWaterPower(10)_:
+0.0001 PowerGeneratingWater(10)
<= 100

c_u_WaterDemand(1)_:
+1 PowerGeneratingWater(1)
<= 9000000

c_u_WaterDemand(2)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
<= 9046340.6499727555

c_u_WaterDemand(3)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
<= 9093721.8804001156

c_u_WaterDemand(4)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
<= 9141727.4591002371

c_u_WaterDemand(5)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
<= 9190149.2699822001

c_u_WaterDemand(6)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
+1 PowerGeneratingWater(6)
<= 9238154.8486823216

c_u_WaterDemand(7)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
+1 PowerGeneratingWater(6)
+1 PowerGeneratingWater(7)
<= 9285119.8469278403

c_u_WaterDemand(8)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
+1 PowerGeneratingWater(6)
+1 PowerGeneratingWater(7)
+1 PowerGeneratingWater(8)
<= 9329327.306968661

c_u_WaterDemand(9)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
+1 PowerGeneratingWater(6)
+1 PowerGeneratingWater(7)
+1 PowerGeneratingWater(8)
+1 PowerGeneratingWater(9)
<= 9372129.9833957627

c_u_WaterDemand(10)_:
+1 PowerGeneratingWater(1)
+1 PowerGeneratingWater(10)
+1 PowerGeneratingWater(2)
+1 PowerGeneratingWater(3)
+1 PowerGeneratingWater(4)
+1 PowerGeneratingWater(5)
+1 PowerGeneratingWater(6)
+1 PowerGeneratingWater(7)
+1 PowerGeneratingWater(8)
+1 PowerGeneratingWater(9)
<= 9414568.4566637538

c_u_MaxStorageCapacity(1)_:
+0 ONE_VAR_CONSTANT
<= 0

c_u_MaxStorageCapacity(2)_:
-1 PowerGeneratingWater(1)
<= -46340.649972755462

c_u_MaxStorageCapacity(3)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
<= -93721.880400115624

c_u_MaxStorageCapacity(4)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
<= -141727.45910023712

c_u_MaxStorageCapacity(5)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
<= -190149.26998220012

c_u_MaxStorageCapacity(6)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
<= -238154.84868232161

c_u_MaxStorageCapacity(7)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
<= -285119.84692784026

c_u_MaxStorageCapacity(8)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
<= -329327.30696866103

c_u_MaxStorageCapacity(9)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
<= -372129.98339576274

c_u_MaxStorageCapacity(10)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
<= -414568.45666375384

c_e_StorageSameAsBefore(1)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(2)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(3)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(4)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(5)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(6)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(7)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(8)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(9)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

c_e_StorageSameAsBefore(10)_:
-1 PowerGeneratingWater(1)
-1 PowerGeneratingWater(2)
-1 PowerGeneratingWater(3)
-1 PowerGeneratingWater(4)
-1 PowerGeneratingWater(5)
-1 PowerGeneratingWater(6)
-1 PowerGeneratingWater(7)
-1 PowerGeneratingWater(8)
-1 PowerGeneratingWater(9)
= -414568.45666375384

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
end
