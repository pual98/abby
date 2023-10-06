#!usr/bin/python3

import math
from calc_coeffs import *

'''
function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
   T   = temperature (degrees C)
   S   = salinity (PSU)
'''


def pHsolver(H):
    '''
    ph = -log([H])
    '''
    return -math.log10(H)

def Hsolver(pH):
    return 10**(-pH)

def OHsolver(H):
    pH=pHsolver(H)
    pOH=14-pH
    OH=10**(-pOH)
    return OH

def CO2solver(k0, pCO2):
    '''
    CO2<->CO2
    k0=[CO2]/pC02
    '''
    return k0*pCO2

def HCO3solver(k1, CO2, H):
    '''
    CO2 + H20 <-> H + HCO3
    k1 = [HCO3][H]/[CO2]
    '''
    return k1*CO2/H

def CO3solver(k2, HCO3, H):
    '''
    HCO3<->H + CO3
    k2 = [CO3][H]/[HCO3]
    '''
    return k2*HCO3/H

def DICsolver(CO2, HCO3, CO3):
    '''
    DIC=[CO2]+[HCO3]+[CO3]
    '''
    DIC=CO2+HCO3+CO3
    return DIC

def BOH4solver(BT,kb, H):
    '''
    [BO4]=BT/(1+[H]/kb)
    '''
    return BT/(1+H/kb)

def TAsolver(HCO3, CO3, BOH4, OH, H):
    '''
    TA = [HCO3] + 2[CO3] + B(OH)4 + [OH] - [H]


    For some reason my answer only matches if I take out these two terms
    +OH
    -H
    '''
    return HCO3+2*CO3+BOH4


def solver(TA, temperature, salinity, pCO2):
    '''
    Write code to solve for pH, [CO2], [HCO3],[CO3], Csat/DIC.

    The following values are known:
    TA (total alkalinity)
    Temperature
    Salinity
    pCO2

    The following values are expected:
    pH
    [CO2]
    [HCO3]
    [CO3]
    Csat

    EX: solver(2350e+6, 20, 34.5, 380e-6)
    '''
    # Use something I don't understand to get me mah ionization constants
    coeffs = calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    BT = coeffs['BT']

    # Get that CO2
    CO2 = CO2solver(k0, pCO2)

    # Random Limits
    pHi=9
    pHlo=6

    TAguess=1

    # Guess ph values until the TA given is close to the calculated TA for a
    # given assumed pH
    while abs(TAguess-TA) > 1e-10:
        H = Hsolver((pHi+pHlo)/2)
        HCO3=HCO3solver(k1, CO2, H)
        CO3=CO3solver(k2, HCO3, H)
        BOH4=BOH4solver(BT, kb, H)
        OH=OHsolver(H)
        TAguess=TAsolver(HCO3, CO3, BOH4, OH, H)
        if TAguess-TA > 0:
            # mean pH is too high, reduce the high spec
            pHi = (pHi+pHlo)/2
        else:
            # mean pH is too low, raise the low spec
            pHlo = (pHi+pHlo)/2

    pH=pHsolver(H)
    DIC = DICsolver(CO2, HCO3, CO3)

    print("ph = " + "{:.5e}".format(pH) + "\n[CO2] = " + "{:.4e}".format(CO2) + "\n[HCO3] = " + \
          "{:.4e}".format(HCO3)+ "\n[CO3] = " + "{:.4e}".format(CO3) + "\nCsat = " + "{:.4e}".format(DIC))

    return

if __name__ == '__main__':
   solver(2350e-6, 20, 34.5, 380e-6)
