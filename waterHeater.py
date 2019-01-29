import math
import ht
import fluids

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def float_range(start, end, step):
    while start <= end:
        yield start
        start += step

###target output is 35C

#Duda Diesel B3-14DW
length1 =  199./1000 #all dims converted to meters
width1 = 86./1000
plateThk1 = 2.35/1000
heatXA1 = .014 #m2
numPlates1 = 20
flowrateAF1 = 5   #GPM

length2 =  199./1000 #all dims converted to meters
width2 = 86./1000
plateThk2 = 2.35/1000
heatXA2 = .014 #m2
numPlates2 = 30
flowrateAF2 = 5   #GPM

#Duda Diesel B3-32DW
length3 =  286./1000 #all dims converted to meters
width3 = 116./1000
plateThk3 = 2.4/1000
heatXA3 = .032 #m2
numPlates3 = 20
flowrateAF3 = 5   #GPM

#constants

tempAF = 65 # Celsius, Webasto temp 78C, minimum
tempWater = 35  # Celsius
# consumeWater = 0.37 #lb/HP/hr
waterTarget = 50

##fluid constants
##AF to 50-50 glycol/water

k_AF = 0.415 # W/(m K), thermal conductivity
C_AF = 3681.9 # J/(kg K), Specific Heat
rho_AF = 1015.6 # kg/m3, Density
mu_AF = 0.000744 # Pa s, Dynamic Viscosity

k_Water = 0.6 # W/(m K), thermal conductivity, water
C_Water = 4184.0 # J/(kg K), Specific Heat, water
rho_Water = 1000. # kg/m3, Density, water
mu_Water = 0.001 # Pa s, Dynamic Viscosity, water

#conversions
# consumeWater /= 2.2 #kg/HP/hr
# consumeWater /= 3600. #kg/HP/s
# consumeWater /= rho_Water #m3/HP/s


def radCalc(numPlates, heatXA, length, width, plateThk, flowrateAF):





    ## calculate surface areas



    areaAF = heatXA * numPlates
    areaWater =  heatXA * numPlates

    areaXSec = plateThk * width
    perimeterXSec = 2. * (plateThk + width)

    ## AF calcs
    flowrateAF = flowrateAF * 3.8 #LPM
    flowrateAF = flowrateAF /60./1000. # convert to m3/s
    massflowAF = flowrateAF *rho_AF
    # print (massflowAF * C_AF)

    ## fluids calcs
    #AF
    Dh_AF = 4. * areaXSec / perimeterXSec
    AFVelocity = flowrateAF/(numPlates / 2.)/(areaXSec)
    reynoldsAF = fluids.core.Reynolds(D=Dh_AF, rho=rho_AF, V=AFVelocity, mu=mu_AF)
    prandltAF = fluids.core.Prandtl(Cp=C_AF , k=k_AF , mu=mu_AF, nu=None, rho=None, alpha=None)
    if reynoldsAF<5000:
        nusseltAF = ht.conv_internal.laminar_Q_const()
    else:
        nusseltAF = ht.conv_internal.turbulent_Dittus_Boelter(Re=reynoldsAF, Pr=prandltAF, heating=False)
    h_AF = nusseltAF * k_AF / Dh_AF


    #setup arrays for data.
    temp = []
    #power = []
    rateWater =[]



    for flowrateWater1 in float_range(2, 4, .1):


        # flowrateWater1 = flowrateWater1 * 3.8 #LPM
        flowrateWater = flowrateWater1 /60./1000. # convert to m3/s
        massflowWater = flowrateWater *rho_Water




        ## fluids calcs
        #Water
        Dh_Water = 4. * areaXSec / perimeterXSec
        WaterVelocity = flowrateWater/(numPlates / 2.)/(areaXSec)
        reynoldsWater = fluids.core.Reynolds(D=Dh_Water, rho=rho_Water, V=WaterVelocity, mu=mu_Water)
        prandltWater = fluids.core.Prandtl(Cp=C_Water , k=k_Water , mu=mu_Water, nu=None, rho=None, alpha=None)
        if reynoldsWater<5000:
            nusseltWater = ht.conv_internal.laminar_Q_const()
        else:
            nusseltWater = ht.conv_internal.turbulent_Dittus_Boelter(Re=reynoldsWater, Pr=prandltWater, heating=True)
        #nusseltWater = 10  # manual overide for laminare flow

        h_Water = nusseltWater * k_Water / Dh_Water



        UA = 1./(h_AF*areaAF) + 1./(h_Water*areaWater)
        UA = 1./UA


        NTU = ht.hx.effectiveness_NTU_method(mc=massflowWater, mh=massflowAF, Cpc=C_Water, Cph=C_AF, subtype='counterflow', Tci=tempWater, Tho=None, Thi=tempAF, Tco=None, UA=UA)
        #NTU = ht.hx.effectiveness_NTU_method(mh=massflowAir, mc=massflowAF, Cph=C_Air, Cpc=C_AF, subtype='crossflow', Thi=tempAir, Tho=None, Tci=tempAF, Tco=None, UA=UA)
        #Power = (NTU['Q']/1000) /0.3/.75
        Tco = (NTU['Tco'])
        Tho = (NTU['Tho'])
        Q = (NTU['Q'])
        #tOut = C_Water / tOut * (rho_Water/1000.) * 15 * 3.8 *10/60.  #sort of calulate time in mins. to heat 15gal water

        temp.append(Tco) #fill data arrays
        rateWater.append((flowrateWater1))

        extra = C_Water * (waterTarget-Tco) * massflowWater
        print (f'flow : {flowrateWater1}  |  main: {Q}  |  watts : {extra}' )

    output = pd.Series(temp, index=rateWater)
    return output



rad1 = radCalc(numPlates1, heatXA1, length1, width1, plateThk1, flowrateAF1)
rad2 = radCalc(numPlates2, heatXA2, length2, width2, plateThk2, flowrateAF2)
# rad3 = radCalc(numPlates3, heatXA3, length3, width3, plateThk3, flowrateAF3)
#rad4 = radCalc(Height4, Width4, Thickness4, fanFlow4, speed4, tempAir4)
#rad5 = radCalc(Height5, Width5, Thickness5, fanFlow5, flowrateAF5, tempAir5)

dataSet = pd.DataFrame({f'{numPlates1} Plate, S'  : rad1,
                        f'{numPlates2} Plate'     : rad2})#,
                        # f'{numPlates3} Plate, L'  : rad3})

                        # f'{numPlates1} Plate'  : rad1,
# print (dataSet)


dataSet.plot()
plt.xlabel('LPM')
plt.ylabel('Celcius')
plt.title('Duda Diesel B3-14DW Comparison')
ymin, ymax = plt.ylim()
plt.text(0.01, ymin + 1, 'Th={0}, Tc={1}'.format(tempAF, tempWater), size=8)
plt.grid(1)
plt.show()

#print dataSet
####testing
