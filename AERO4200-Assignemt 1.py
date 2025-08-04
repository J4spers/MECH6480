#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:25:21 2024
    
@author: Jasper Macqueen Wetzel
"""
import math
import numpy
from pylab import plot,show,legend, xlabel, ylabel, grid, xlim, ylim, title
from sympy.solvers import solve
from sympy import Symbol
v = Symbol("v")
from scipy.optimize import fmin
from Air_value_finder import density
import pytexit


## Assumptions ##
g = 9.81
Density_Cruise = density(2500) ## Air Denisty at 2500m [kg/m^3] ##
Temp_Cruise = 271.92 ## Temparcture at 2500, [K]elvin ##
Pressure_Cruise_inf = 7.4692 ## Surounding Air Presure at cruise [N/m^2] ##
Viscosity_Cruise = 1.70995064 * 10**(-5) ## ? check this ##

## List of given Parameters ##

Wingspan = 4 ##[m], b ##
Cruise_Wing_Loading = 40 ##[N/m^2]##
Max_Load_Factor = 2.1 
V_Cruise = 60/3.6 ## divided by 3.6 to get to m/s from km/h, [m/s]##
Alt_Cruise = 2500 ##[m]##
e = 0.94 ## Span and Oswald efficency factor ##
Mass = 12.5 ## [kg] ##
Payload = 2 ## [kg] ##
Weight = (Mass + Payload)*g ## include the payload in the weight to acount for the payload in design ##

Engine_Power = 450 *2 ## [W]atts  Note: there are 2 engines, each engine has 450 W of power## 
Engine_Eff_sea = 0.84  ## Efficency reduces Linearly ##
Engine_Eff_6000 = 0.42 ## at 6000m altitude ##

Cd_fuselarge = 0.024 ## uses wing plan area as lenght scale ##

## AeroFoil Data NACA 2412 ##

Cd_0_2412 = 0.006 ## taken from tables as drag at zero lift at cruise Re##

"""
    Airofoil tools gives Cd0 at Re 1,000,000 as 0.0568
    Airofil tools gives Cd0 at Re 500,000
"""

def finite_Lift_Slope(a0, e, Wing_AR):
    a = (a0)/(1+(a0/(math.pi*e*Wing_AR)))
    return a

def Re(Density, Velocity, Length, Dynamic_Viscosity):
    Reynolds_number = (Density * Velocity * Length)/Dynamic_Viscosity
    return Reynolds_number



Step_Size = 120*100


print("""
Question 1. 
Cruise parameters

""")

L_Cruise = Weight ## In trim the cruise lift is equal to the weight of the aircraft ##

Wing_Area = Weight/Cruise_Wing_Loading ## W/S, where S is the Wing Area, [m] ##
Wing_AR = (Wingspan**2)/Wing_Area ## Wing aspect ratio (AR), b^2 /S ##
Cord_L = Wing_Area/Wingspan

a_0 = ((-1.0)-(1.2))/(10 - (-10)) *-1
a = finite_Lift_Slope(a_0,e,Wing_AR)

Cl_Cruise = L_Cruise/(Wing_Area * 0.5 * Density_Cruise * V_Cruise**2)

AoA_trim = -2 ## [Deg] for a Cl of 0.024 the NACA 2412 table gives a AoA of aproximetly -2 deg ##

Abs_AoA = abs(AoA_trim) + abs(Cl_Cruise/a)

Cd_0 = Cd_fuselarge + Cd_0_2412 ## combine to find Cd_0 of aircaft ##

CD_Cruise = Cd_0 + (Cl_Cruise**2 /(math.pi*e*Wing_AR))

D_Cruise = CD_Cruise * Wing_Area * 0.5 * Density_Cruise * V_Cruise**2 ## [N] ##

Re_Cruise = Re(Density_Cruise,V_Cruise,Cord_L,Viscosity_Cruise)

Cl_Max = 1.5 ## From Tables @ aprox 15 deg ##
"""
** maybe remove this ***
    At a Renolds Number of 200,000 Cl max was recorded to be 1.3192 at 14.25 deg
    At a Renolds Number of 500,000 Cl max was recorded to be 1.4070 at 15.25 deg
    At a Renolds Number of 1,000,000 Cl ma was recorded to be 1.5820 at 16.5 deg
    Given a Cruise Renolds number of 
    
    NACA data gives Cl max as aprox 1.5
"""


def V_Stall(Altitude,Weight,Wing_Area):
    V_Stall = ((2*Weight)/(density(Altitude) * Wing_Area * Cl_Max))**0.5
    return V_Stall

def Cl(Lift, Altitude, Velocity):
    return Lift/(Wing_Area * 0.5 * density(Altitude) * Velocity**2)

def CD(Cl,Wing_AR):
    return Cd_0 + (Cl**2 /(math.pi*e*Wing_AR))

def Lift(Cl, Altitude, Velocity):
    return Cl * Wing_Area * 0.5 * density(Altitude) * Velocity**2
def Drag(CD,Altitude, Velocity):
    return CD * Wing_Area * 0.5 * density(Altitude) * Velocity**2

print("Lift:",Lift(Cl_Cruise, Alt_Cruise, V_Cruise), '\n'  )
print("Drag:", Drag(CD_Cruise,Alt_Cruise,V_Cruise),'\n')
print("Absolute Angle of Attack:",Abs_AoA,'\n')


print("""
Question 2.
Power Required vs Power Avalible

""")


def Engine_eff(Altitude):
    Engine_eff = ((Engine_Eff_sea - Engine_Eff_6000)/(0-6000)) * Altitude + Engine_Eff_sea
    return Engine_eff

def Power_Required(Velocity, Altitude, Wing_Area, Wing_AR): #add altitude ie make independent of cruise peramiters ##
    Power_required = (0.5* density(Altitude) * (Velocity**3) * Wing_Area*Cd_0) + (((Weight**2/(0.5*density(Altitude)*Velocity*Wing_Area)))/(math.pi * e *Wing_AR))
    return Power_required

def Power_Avalible(Altitude):  
    Power_Avalible = Engine_eff(Altitude) * Engine_Power
    return Power_Avalible ## returns units [W] ##

Velocity_range = numpy.linspace(0.01,120,Step_Size) ## from 0m/s to speed of sound at sea level, under the assumption the plane will travel at sub-sonic speeds ##
Altitude_Range = numpy.linspace(1,12000,Step_Size)
"""
    The Velocity Range was initial set from 1 m/s to the speed of sound 
    at sea level (343 m/s). Once the glider was found to not exceed aprox. 30 m/s
    the Velocity range was decresed to 1-120 m/s and the step size increased to 100
    to increase precision and remove unessiary computaions.
"""

Required_Cruise_Power = Power_Required(Velocity_range, Alt_Cruise,Wing_Area,Wing_AR)

plot(Velocity_range,Required_Cruise_Power, label = 'Required Power')
plot(Velocity_range, Power_Avalible(2500)*Velocity_range*(1/Velocity_range),'r-' ,label = 'Power Available')
legend()

xlabel("Velocity [m/s]"),ylabel("Power [W]"),grid(),title("Required Cruise Power at 2500m")

xlim(0,30)
ylim(0,2*10**3)

## power avalible is suitabel ##
print("Power Required at Cruise:",Power_Required(V_Cruise, Alt_Cruise,Wing_Area,Wing_AR),'\n')
print("Power Avalible at Cruise:", Power_Avalible(Alt_Cruise),'\n')

print("""
Question 3.
Maximum Speed 

""")

Max_Speed = solve(Engine_eff(2500)*Engine_Power - (0.5*Density_Cruise * (v**3) * Wing_Area*Cd_0) + (Weight**2/(0.5*Density_Cruise*v*Wing_Area))/(math.pi * e *Wing_AR)) ## finds the intersect of required power and power avalible ## 
print("The maximum possible speed at 2500m is",Max_Speed[1],"[m/s]")

Cruise_Alt_Velocity_Range = numpy.linspace(7.5,25,35)

print("""
Question 4.
Endurance Possible at Cruising Conditions 
    and
Maximum Endurance Possible at the Cruising Altitude (2500m)
""")

Battery_charge = 6.7 * 4 ## [A/h]  multiplied by the number of batteries##
Voltage = 14.8

def Endurance(Altitude, Velocity,Wing_Area,Wing_AR): ## [h] ##
    Flight_time = [(Battery_charge*Voltage*Engine_eff(Altitude))/Power_Required(Velocity,Altitude,Wing_Area,Wing_AR),Velocity]
    return Flight_time

"""
     a) 
"""
Cruise_Endurance = Endurance(Alt_Cruise, V_Cruise,Wing_Area,Wing_AR)[0] ## [Hours] ##

"""
     b) 
"""

Max_Cruise_Endurance = Endurance(Alt_Cruise, Velocity_range,Wing_Area,Wing_AR) ## [Hours] ##


print("The Cruise endurance is",Cruise_Endurance,'\n')
print("The maximum endurance is",max(Max_Cruise_Endurance[0]), "Hours, with a velocity of",
      Max_Cruise_Endurance[1][numpy.where(Max_Cruise_Endurance[0]==max(Max_Cruise_Endurance[0]))[0]],"m/s")


print("""
Question 5,
Possible Range at Cruise 

""")
def Range(Altitude, Velocity,Wing_Area,Wing_AR):
    Range = Velocity * Endurance(Altitude, Velocity,Wing_Area,Wing_AR)[0] * 3600
    return [Range, Velocity]

"""
     a) 
"""
Cruise_Range = Range(Alt_Cruise,V_Cruise,Wing_Area,Wing_AR) ## [m] ##

"""
     b) 
"""
Max_Cruise_Range = Range(Alt_Cruise,Velocity_range,Wing_Area,Wing_AR) ## [m] ##\

print("The Range at Cruise Condions is:",Cruise_Range,'\n')
    
print("The maximum range is",(max(Max_Cruise_Range[0]))/1000, "km, at a velocity of",
      Max_Cruise_Range[1][numpy.where(Max_Cruise_Range[0]==max(Max_Cruise_Range[0]))[0]],"m/s")
print("""
Question 6,
Absiolute Ceiling,
saved value: 4243.0 [m]
Tolerance = 10
""")
def Absolute_Ceiling_Finder(Tolerance):
    Absolute_Ceiling = []
    for i in Altitude_Range:
        if abs(Power_Avalible(i) - min(Power_Required(Velocity_range,i))) < Tolerance:
            Absolute_Ceiling.append(i)
        else:
            print(i)
            continue
        
    return max(Absolute_Ceiling)

Absolute_Ceiling = 8928

"""
    This function finds the diffrence between the Power Required and the Power Avalible
    for velocities from 1 m/s to 120 m/s at diffrent altitudes. Starting at 1m and continuing 
                    till 6000m, calculating the diffrence ever 0.5m. 
    It then saves any values within a given tolerance, returning the highst value.
"""
print("The absolute Ceiling of the Aircarft is" , Absolute_Ceiling, "meters",'\n')

print("""
Question 7,
Maximum Level Turn Rate and Minimum Radius of Curviture 

""")            

def Turn_Rate(Load_Factor,Altitude): ## [rad/s] ##
    return g*((0.5 * (g**2)*Cl_Max*density(Altitude)*Wing_Area*((Load_Factor**2)-1))/(Weight*Load_Factor))**(0.5) #(g*(((Load_Factor**2)-1)**(0.5)))/(Velocity)

def Turn_Rate_Velocity(Load_Factor, Velocity):
    return g*(((Load_Factor**2) -1)**(0.5))/Velocity

def Turn_Velocity(Turn_Rate, Load_Factor):
    return g*(((Load_Factor**2) -1)**(0.5))/Turn_Rate
      
def Turn_Radius(Load_Factor,Velocity):
    return ((Velocity**2))/(g*(((Load_Factor**2)-1)**(0.5)))

def Min_Power_Speed(Altitude):
    Min_Power_Speed = solve(Engine_eff(Altitude)*Engine_Power - (0.5* density(Altitude) * (v**3) * Wing_Area*Cd_0) + (Weight**2/(0.5* density(Altitude                                                                                                                                   ) *v*Wing_Area))/(math.pi * e *Wing_AR))
    return abs(Min_Power_Speed[0])


"""
    To Maximise turn radius minimise velocity at each condition,
            however velocity is set for cruise conditions.
"""

def Level_Turn(Load_Factor, Altitude, Weight, Cl):
    Turn = Turn_Rate(Load_Factor, Altitude)
    
    if Turn_Velocity(Turn, Load_Factor) < V_Stall(Altitude, Weight,Wing_Area):
        Level_Turn = [(g/(1.2 * V_Stall(Altitude, Weight,Wing_Area)))  ,   
                      ((1.2*V_Stall(Altitude, Weight,Wing_Area))**2)/g]
    else:
        Level_Turn = [Turn,Turn_Radius(Load_Factor, Turn_Velocity(Turn, Load_Factor))]
        
    return Level_Turn

Turn_Max_Cruise = Level_Turn(Max_Load_Factor,Alt_Cruise,Weight,Cl_Max)
Turn_Max_Sea = Level_Turn(Max_Load_Factor, 0, Weight, Cl_Max)
Turn_Max_Ceiling = Level_Turn(Max_Load_Factor, Absolute_Ceiling, Weight, Cl_Max)

print("The maximum turn rate at cruise altitude is", Turn_Max_Cruise[0],"rad/sec, and has an assoicated radius of",Turn_Max_Cruise[1],"[m]",'\n')
print("The maximum turn rate at Sea altitude is", Turn_Max_Sea[0],"rad/sec, and has an assoicated radius of",Turn_Max_Sea[1],"[m]",'\n')
print("The maximum turn rate at Ceiling altitude is", Turn_Max_Ceiling[0],"rad/sec, and has an assoicated radius of",Turn_Max_Ceiling[1],"[m]",'\n')


"""
    If the Turn Velocity is found to be bellow the stall velocity the actual minimum turning velocity is at V-Stall.
    Using circular motion mechnaics as shown in "MECH3410 Section C2 Tutorial", V^2 /rg = tan( Theta ) can be used 
            to find the radius (r) by optimising for theta. Then omega = v/r can be substituted to get,
                                                w = g/v and r = v^2/g.
                        Since turning will stall one of the wings, use v =1.2 * V Stall.
    
"""



print("""
Question 8,
""")
## Provided Values ##
ICE_Eff_sea = 0.8

ICE_Eff_5000 = 0.4
ICE_Power = 480 * 2 ## [W] ##
ICE_Fuel_Rate = 0.63*2/(1000*60) ## [kg/W/min] ##
Density_of_Kerosine = 0.82*1000## [kg/m^3] ##
Fuel_Capacity = 12 /1000 ## [m^3] ##
ICE_mass = 13.6 ## [kg] ##
ICE_Weight_Cruise = (ICE_mass + Payload + 0.6 * Fuel_Capacity * Density_of_Kerosine)*g
ICE_Weight_Max = (ICE_mass + Payload + (Fuel_Capacity*Density_of_Kerosine)) * g
ICE_Weight_Min = (ICE_mass + Payload) * g
Cl_Cruise_Average = (Cl(ICE_Weight_Max,2500,V_Cruise) - Cl(ICE_Weight_Min,2500,V_Cruise))/2
## assume Wing parametrs do not change ##
ICE_Wing_Area = ICE_Weight_Cruise/40 ## [m^2] ##
ICE_Wing_AR = (Wingspan**2)/ICE_Wing_Area

def ICE_eff(Altitude):
    ICE_eff = ((ICE_Eff_sea - ICE_Eff_5000)/(0-5000)) * Altitude + ICE_Eff_sea
    return ICE_eff

def ICE_Power_Required(Velocity, Altitude, Weight, Wing_Area, Wing_AR): #add altitude ie make independent of cruise peramiters ##
    ICE_Power_required = (0.5* density(Altitude) * (Velocity**3) * Wing_Area*Cd_0) + (((Weight**2/(0.5*density(Altitude)*Velocity*Wing_Area)))/(math.pi * e *Wing_AR))
    return ICE_Power_required

def ICE_Power_Avalible(Altitude):  
    ICE_Power_Avalible = ICE_eff(Altitude) * ICE_Power
    return ICE_Power_Avalible ## returns units [W] ##

def Inst_Fuel_Consumption(Velocity,Altitude,Wing_Area,Wing_AR): ## returns flight time in Hrs ##
    Minutes = 0
    w0 = ICE_Weight_Max
    while w0 > (ICE_Weight_Min):
        dw = ICE_Fuel_Rate * (ICE_Power_Required(Velocity, Altitude, w0, Wing_Area, Wing_AR)/ICE_eff(Altitude)) * g
        Minutes = Minutes +1
        w0 = w0 - dw
        print(w0)
    return Minutes/60

"""
    Returns Endurnace at Cruise conditions as 10.23 hrs
"""

def ICE_Endurance_Max(Altitude, Weight, Wing_AR ):
    ICE_Endurance_Max = (ICE_eff(Altitude)/(ICE_Fuel_Rate * g)) *   (((3*Cd_0* math.pi * e * Wing_AR)**(3/4))/(4*Cd_0))  * ((density(Altitude))*2)**(0.5) * (ICE_Weight_Min**(-0.5) - Weight**(-0.5))
    return ICE_Endurance_Max/60
    

def ICE_Range(Altitude,Velocity,Weight,Wing_Area,Wing_AR):
    Cl = ICE_Weight_Cruise/(0.5*density(Altitude)*Velocity**2 * Wing_Area)
    Cd =CD(Cl,ICE_Wing_AR)
    Range = (ICE_eff(Altitude)/(ICE_Fuel_Rate * g/60))* (Cl/Cd) * numpy.log(Weight/ICE_Weight_Min)
    
    return Range/1000

def ICE_Range_Max(Altitude,Weight,Wing_AR):
    Max_Range = (ICE_eff(Altitude)/(ICE_Fuel_Rate * g/60))* ((math.pi * e * Wing_AR)/(4*Cd_0))**(0.5) * numpy.log(Weight/ICE_Weight_Min)
    return Max_Range/1000

ICE_Required_Cruise_Power = ICE_Power_Required(V_Cruise, Alt_Cruise, ICE_Weight_Cruise,ICE_Wing_Area,ICE_Wing_AR)
ICE_Endurnace = 10.23
ICE_Range_Cruise = (ICE_Endurnace *60 *60 * V_Cruise)/1000 ## [km] ##

print("ICE at Cruise",'\n',
      "Cruise Endurance", ICE_Endurnace,'\n',
      "Cruise Range", ICE_Range_Cruise,'\n')

print("""
Question 9,
""")
New_WingSpan = Wingspan * 2     
New_WingArea_ICE = ICE_Weight_Cruise/62 ## [m^2] ##
New_WingArea = Weight/62 ## [m^2] ##
New_Wing_AR = (New_WingSpan**2)/New_WingArea
New_Wing_AR_ICE =(New_WingSpan**2)/New_WingArea_ICE

"""
    Internal Combustion Engine
"""
#New_ICE_Endurance = Inst_Fuel_Consumption(V_Cruise,Alt_Cruise,New_WingArea_ICE,New_Wing_AR_ICE) ## [Hours] ##
    ## ^ only run if required, returns 18.5 Hours ##
ICE_E =18.5

New_ICE_Range = ICE_Range(Alt_Cruise,V_Cruise,ICE_Weight_Max,New_WingArea_ICE,New_Wing_AR_ICE) ## [km] ##

New_ICE_MaxEndurance = ICE_Endurance_Max(Alt_Cruise, ICE_Weight_Max, New_Wing_AR_ICE) ## [Hours] ##

New_ICE_MaxRange = ICE_Range_Max(Alt_Cruise,ICE_Weight_Max,New_Wing_AR_ICE) ## [km] ##

print("New Wing, ICE at Cruise",'\n',
      "Cruise Endurance", ICE_E,'\n',
      "Cruise Range", New_ICE_Range,'\n')


"""
    Electric Motor
"""
New_Cruise_Endurance = Endurance(Alt_Cruise, V_Cruise, New_WingArea, New_Wing_AR)[0] ## [h] ##

New_Cruise_Range = (Range(Alt_Cruise,V_Cruise,New_WingArea,New_Wing_AR)[0])/1000 ## [km] ##

New_MaxEndurance = max(Endurance(Alt_Cruise, Velocity_range, New_WingArea, New_Wing_AR)[0]) ## [h] ##

New_MaxRange = max(Range(Alt_Cruise,Velocity_range,New_WingArea,New_Wing_AR)[0])/1000 ## [km] ##

print("New Wing, Electric Motor at Cruise",'\n',
      "Cruise Endurance", New_Cruise_Endurance,'\n',
      "Cruise Range", New_Cruise_Range,'\n')


print("""
Question 10,
""")

print("Original Wing, Electric Motor at Cruise",'\n',
      "Cruise Endurance", Cruise_Endurance,'\n',
      "Cruise Range", Cruise_Range[0]/1000,'\n')

print("Recomended Improvment, Electric Motor at Cruise",'\n',
      "Cruise Endurance", Endurance(Alt_Cruise, 50/3.6, Wing_Area, Wing_AR)[0],'\n',
      "Cruise Range", (Range(Alt_Cruise,50/3.6,Wing_Area,Wing_AR)[0])/1000,'\n')

print('\n',"As can be seen by the Values Tabulated above, the New Wing with an", '\n' ,
      "Internal Combustion Engine provides the largest increase to Endurance and Range", '\n')


print("""
Questiomn 11,
""")
## NACA 65-210 ##

CM_0 = 0.07
#Cg_Pos = 0.4 * Cord_L
#Static_Margin = 0.15 * Cord_L
Cg_Pos = 0.4
Static_Margin = 0.15
l_t = 1.4
e_tail = 0.86
np_Pos = Cg_Pos + Static_Margin
d3_dalpha = 0.32
ac = 0.25 * Cord_L ## from NACA tables

#Lochie:
#a_0 = (0.6-(-0.5))/(4-(-8))
#a_t_0 = (0.4-(-0.6))/(4-(-6))

a_T_0 = ((-0.8)-(0.8))/(8 - (-8)) * -1

SC = 3.5

a_T = finite_Lift_Slope(a_T_0, e_tail, SC )


V_H = ((np_Pos-ac)) / ((a_T/a)*(1   -  d3_dalpha))
S_t = (SC * V_H)/l_t

CM_ac_wb = -0.05
Tail_Setting_Angle = (((CM_0-CM_ac_wb)/(V_H*a_T))) ## [Deg] ##

# Jasper:
# """
#     Lift Slopes of Tail and Wing collected from NACA tables
# """
# a_T_0 = ((-0.8)-(0.8))/((-8*(math.pi/180)) -   (8*(math.pi/180)))
# a_0 = ((-1.0)-(1.2))/((-10*(math.pi/180)) -   (10*(math.pi/180)))

#     ## Characteristic scale should be between 2 and 5, let SC = 2.25 ##
# SC = 3.5

# def finite_Lift_Slope(a0, e, Wing_AR):
#     a = (a0)/(1+(a0/(math.pi*e*Wing_AR)))
#     return a

# a_T = finite_Lift_Slope(a_T_0, e_tail, SC )
# a = finite_Lift_Slope(a_0,e,Wing_AR)

# V_H = ((np_Pos - ac)) / ((a_T/a)*(1   -  d3_dalpha))


# S_t = (SC * V_H)/l_t

# """
#       Aerodynamics of the wing-body are dominated by those of the wing 
#                   so CM_AC = CM_AC of the wing body.
#                 From Tables this gives CM_AC = -0.5
# """
# CM_ac_wb = -0.05
# def Downwash(angle):
#     return 0.0+ 0.32 * (angle *(math.pi/180))

# Tail_Setting_Angle = (((CM_0-CM_ac_wb)/(V_H*a_T)) - Downwash(0)) ## [Deg] ##

print("Tail Plane Area:",S_t,'\n')
print("Tail Plane Volume:",V_H,'\n')
print("setting angle",Tail_Setting_Angle,'\n')


print("""
Question 12,
""")

"""
        Static Stability Check
"""
Alpha_e = -0.07 /(a*(    (((Cg_Pos) - (ac))) -  (V_H* (a_T/a)* (1-d3_dalpha)   )))  
print("Alpha_e",Alpha_e,'\n')
print("dCM/dAlpha", a*(    ((Cg_Pos - ac)) -  V_H* (a_T/a)* (1-d3_dalpha)   ),'\n')
print("Check for negative moment slope:" ,0  >   a*(    ((Cg_Pos - ac)) -  V_H* (a_T/a)* (1-d3_dalpha)   ),'\n'  ,'\n'                         )
print("Check for possitive moment at zero angle of attach :" , 0 <      CM_ac_wb  +  V_H*a_T*Tail_Setting_Angle,'\n'  ,'\n' )
print("Check alpha  < stall angle :", Abs_AoA < 16.5,'\n'  ,'\n' )

print("""
Question 13
""")
Elevator_Control_Eff = 0.04 ## [per deg] ##

Trim = (CM_0 + (a*(    ((Cg_Pos - ac)) -  V_H* (a_T/a)* (1-d3_dalpha)   )*(Abs_AoA)*(math.pi/180)))/(V_H*Elevator_Control_Eff)

print("The tail trim angle at cruise is:",Trim,"degrees",'\n')

print("""
Question 14
Runway length

""")
def V_Landing(Altitude): 
    V_Landing = 1.3 * V_Stall(Altitude, Weight, Wing_Area)
    return V_Landing ## [m/s] ##

Cl_Max_Landing = 1.4
Friction_Factor = 0.13

"""
    Assume the change in air densty during landing is negligable
    Assume Temprature effects are negligable
"""

def S_Landing(Altitude,Weight):
    Lift = 0.5*density(Altitude)* ((V_Landing(Altitude)*0.7)**2) * Wing_Area * Cl_Max_Landing
    Drag = 0.5*density(Altitude)* ((V_Landing(Altitude)*0.7)**2) * Wing_Area * Cd_0
    
    Runway_Length = -(1.69 * Weight**2)/(density(Altitude)*(-g)*Wing_Area*Cl_Max_Landing*(Drag + (Friction_Factor*(Weight))))
    return Runway_Length ## [m] ##

Landing_Length_Sea_Level = S_Landing(0,Weight)
Landing_Length_1500m = S_Landing(1500,Weight)

print("The runway length at sea level is",Landing_Length_Sea_Level,'\n')
print("The runway length at 1500m is",Landing_Length_1500m,'\n')

print("""
Question 15
""")      
#### !!!! ###
