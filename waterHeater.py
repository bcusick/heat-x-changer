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
numPlates1 = 10
flowrateCoolant1 = 5.   #GPM

length2 =  199./1000 #all dims converted to meters
width2 = 86./1000
plateThk2 = 2.35/1000
heatXA2 = .014 #m2
numPlates2 = 30
flowrateCoolant2 = 5.   #GPM

#Duda Diesel B3-32DW
length3 =  286./1000 #all dims converted to meters
width3 = 116./1000
plateThk3 = 2.4/1000
heatXA3 = .032 #m2
numPlates3 = 20
flowrateCoolant3 = 5.   #GPM


#constants

tempCoolant = -5 # Celsius, Webasto temp 78C, minimum
tempFuel = 30  # Celsius
consumeFuel = 0.37 #lb/HP/hr


##fluid constants
##coolant to 50-50 glycol/water

k_Coolant = 0.415 # W/(m K), thermal conductivity
C_Coolant = 3681.9 # J/(kg K), Specific Heat
rho_Coolant = 1015.6 # kg/m3, Density
mu_Coolant = 0.000744 # Pa s, Dynamic Viscosity

k_Fuel = 0.6 # W/(m K), thermal conductivity, water
C_Fuel = 4184.0 # J/(kg K), Specific Heat, water
rho_Fuel = 1000. # kg/m3, Density, water
mu_Fuel = 0.001 # Pa s, Dynamic Viscosity, water

#conversions
consumeFuel /= 2.2 #kg/HP/hr
consumeFuel /= 3600. #kg/HP/s
consumeFuel /= rho_Fuel #m3/HP/s


def radCalc(numPlates, heatXA, length, width, plateThk, flowrateCoolant):





    ## calculate surface areas



    areaCoolant = heatXA * numPlates
    areaFuel =  heatXA * numPlates

    areaXSec = plateThk * width
    perimeterXSec = 2. * (plateThk + width)

    ## coolant calcs
    flowrateCoolant = flowrateCoolant * 3.8 #LPM
    flowrateCoolant = flowrateCoolant /60./1000. # convert to m3/s
    massflowCoolant = flowrateCoolant *rho_Coolant
    print massflowCoolant * C_Coolant

    ## fluids calcs
    #Coolant
    Dh_Coolant = 4. * areaXSec / perimeterXSec
    coolantVelocity = flowrateCoolant/(numPlates / 2.)/(areaXSec)
    reynoldsCoolant = fluids.core.Reynolds(D=Dh_Coolant, rho=rho_Coolant, V=coolantVelocity, mu=mu_Coolant)
    prandltCoolant = fluids.core.Prandtl(Cp=C_Coolant , k=k_Coolant , mu=mu_Coolant, nu=None, rho=None, alpha=None)
    if reynoldsCoolant<5000:
        nusseltCoolant = ht.conv_internal.laminar_Q_const()
    else:
        nusseltCoolant = ht.conv_internal.turbulent_Dittus_Boelter(Re=reynoldsCoolant, Pr=prandltCoolant, heating=False)
    h_Coolant = nusseltCoolant * k_Coolant / Dh_Coolant


    #setup arrays for data.
    temp = []
    #power = []
    rateFuel =[]



    for flowrateFuel1 in float_range(.01, 0.5, .01):


        flowrateFuel = flowrateFuel1 * 3.8 #LPM
        flowrateFuel = flowrateFuel /60./1000. # convert to m3/s
        massflowFuel = flowrateFuel *rho_Fuel




        ## fluids calcs
        #Fuel
        Dh_Fuel = 4. * areaXSec / perimeterXSec
        fuelVelocity = flowrateFuel/(numPlates / 2.)/(areaXSec)
        reynoldsFuel = fluids.core.Reynolds(D=Dh_Fuel, rho=rho_Fuel, V=fuelVelocity, mu=mu_Fuel)
        prandltFuel = fluids.core.Prandtl(Cp=C_Fuel , k=k_Fuel , mu=mu_Fuel, nu=None, rho=None, alpha=None)
        if reynoldsFuel<5000:
            nusseltFuel = ht.conv_internal.laminar_Q_const()
        else:
            nusseltFuel = ht.conv_internal.turbulent_Dittus_Boelter(Re=reynoldsFuel, Pr=prandltFuel, heating=True)
        #nusseltFuel = 10  # manual overide for laminare flow

        h_Fuel = nusseltFuel * k_Fuel / Dh_Fuel



        UA = 1./(h_Coolant*areaCoolant) + 1./(h_Fuel*areaFuel)
        UA = 1./UA


        NTU = ht.hx.effectiveness_NTU_method(mh=massflowFuel, mc=massflowCoolant, Cph=C_Fuel, Cpc=C_Coolant, subtype='counterflow', Thi=tempFuel, Tho=None, Tci=tempCoolant, Tco=None, UA=UA)
        #NTU = ht.hx.effectiveness_NTU_method(mh=massflowAir, mc=massflowCoolant, Cph=C_Air, Cpc=C_Coolant, subtype='crossflow', Thi=tempAir, Tho=None, Tci=tempCoolant, Tco=None, UA=UA)
        #Power = (NTU['Q']/1000) /0.3/.75
        tOut = (NTU['Q'])
        #tOut = C_Fuel / tOut * (rho_Fuel/1000.) * 15 * 3.8 *10/60.  #sort of calulate time in mins. to heat 15gal water

        temp.append(tOut) #fill data arrays
        rateFuel.append(flowrateFuel1)

    output = pd.Series(temp, index=rateFuel)
    return output



rad1 = radCalc(numPlates1, heatXA1, length1, width1, plateThk1, flowrateCoolant1)
rad2 = radCalc(numPlates2, heatXA2, length2, width2, plateThk2, flowrateCoolant2)
rad3 = radCalc(numPlates3, heatXA3, length3, width3, plateThk3, flowrateCoolant3)
#rad4 = radCalc(Height4, Width4, Thickness4, fanFlow4, speed4, tempAir4)
#rad5 = radCalc(Height5, Width5, Thickness5, fanFlow5, flowrateCoolant5, tempAir5)

dataSet = pd.DataFrame({'{0} Plate'.format(numPlates1)      : rad1,
                        '{0} Plate'.format(numPlates2)      : rad2,
                        '{0} Plate'.format(numPlates3)      : rad3})
print dataSet

dataSet.plot()
plt.xlabel('GPM')
plt.ylabel('Celcius')
plt.title('Duda Diesel B3-14DW Comparison')
ymin, ymax = plt.ylim()
plt.text(0.01, ymin + 1, 'Th={0}, Tc={1}'.format(tempCoolant, tempFuel), size=8)
plt.grid(1)
plt.show()

#print dataSet
####testing
