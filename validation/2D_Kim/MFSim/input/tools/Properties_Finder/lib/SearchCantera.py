#!/usr/bin/env python
# -*- coding:utf-8 -*-

#'''
#--------------------------------------
#Fluid Properties Finder Library
#
#Global Variables:
#`user_param`
#`fluid`
#
#Basic Routines:
#`start`
#`end`
#`input`
#`set_cantera`
#`get_cantera`
#
#--------------------------------------
#'''

# Import system libraries
import commentjson
import math
import cantera as ct

# Global variables
user_param = { # Default values
    "Substance":{
        "Composition":"O2:1.0 , N2:3.76",
        "Fraction Type":"molar"
    },
    "Termodinamic State":{
        "Temperature": 20.0, # [°C]
        "Pressure": 101325.0 # [Pa]
    },
    "Properties":[
        "rho", # [kg/m³]
        "mu", # [Pa s]
        "cp", # [J/(kg K)]
        "k" # [W/(m K)]
    ],
    "cti_database":"lib/gri30.cti"
}
fluid = 0

#::::::::::::::::::::::::::::::::::::::::::::::::::
def start():
    print("------------------------------------")
    print("PropFind - A Cantera API")
    print("------------------------------------")
    
#::::::::::::::::::::::::::::::::::::::::::::::::::
def end():
    print("")
    print("------------------------------------")
    print("Finished")
    print("------------------------------------")

#::::::::::::::::::::::::::::::::::::::::::::::::::
def input():
    global user_param
    print('Reading input file ...')
    
    with open('Fluid_input.jsonc') as jsonc_file:
        json_input = commentjson.load(jsonc_file)

    user_param.update(json_input)

#::::::::::::::::::::::::::::::::::::::::::::::::::
def set_cantera():
    global user_param
    global fluid
    print('Setting info into Cantera library ...')
    
    #Entry of cti file
    fluid = ct.Solution( user_param["cti_database"] )

    #Entry of substance composition
    if (user_param["Substance"]["Fraction Type"]=="molar" ):
        fluid.X = user_param["Substance"]["Composition"]
    elif(user_param["Substance"]["Fraction Type"]=="mass" ):
        fluid.Y = user_param["Substance"]["Composition"]
    
    Temp = user_param["Termodinamic State"]["Temperature"] + 273.15 # Transform it to K
    Pres = user_param["Termodinamic State"]["Pressure"]
    
    #Enty of Pressure and Temperature
    fluid.TP = Temp, Pres

#::::::::::::::::::::::::::::::::::::::::::::::::::
def get_cantera():
    global fluid
    print('Getting Properties ...')
    print('')
    
    print('Printing Properties for:')
    print('Fluid: "{0:s}"'.format(user_param["Substance"]["Composition"]))
    print('In the Temperature of {0:f} °C'.format(user_param["Termodinamic State"]["Temperature"]))
    print('And Pressure of {0:f} Pa'.format(user_param["Termodinamic State"]["Pressure"]))
    print('')

    if("rho" in user_param["Properties"]):
        print( '{0:25s} | {1:12E} | kg/m³'.format('Density',fluid.density) )

    if("mu" in user_param["Properties"]):
        print( '{0:25s} | {1:12E} | Pa s'.format('Dynamic Viscosity',fluid.viscosity) )

    if("cp" in user_param["Properties"]):
        print( '{0:25s} | {1:12E} | kJ/(kg K)'.format('Specific Heat',fluid.cp_mass) )

    if("k" in user_param["Properties"]):
        print( '{0:25s} | {1:12E} | W/(m K)'.format('Thermal Conductivity',fluid.thermal_conductivity) )

    if("c" in user_param["Properties"]):
        gamma = fluid.cp_mass/fluid.cv_mass
        R = ct.gas_constant/fluid.mean_molecular_weight
        c = math.sqrt(gamma*R*fluid.T)
        print( '{0:25s} | {1:12E} | m/s'.format('Sound Speed',c) )

#::::::::::::::::::::::::::::::::::::::::::::::::::
