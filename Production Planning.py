# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:11 2021

@author: user
"""
import numpy as np
import math


def mad(estimates,actual_data):
    mad = np.mean(abs(estimates-actual_data))
    return mad


def moving_averages(int_data,target_data,n):
    estimates= []
    for i in range(len(target_data+1)):
        if i<n:
            estimate = np.mean(np.append(int_data[-n+i:],target_data[0:i]))
        else:
            estimate = np.mean(target_data[i-n:i])
        estimates.append(estimate)
    estimates = np.array(estimates)
    return estimates

def seasonal_stationary(int_data,n, season_length):
    avg =np.mean(int_data)
    seasonal_factors = int_data / avg
    avg_factors = np.zeros(season_length)
    seasons = 0
    for i in range (len(int_data)//season_length):
        avg_factors += seasonal_factors[season_length * i: season_length*i+season_length] 
        seasons +=1
    avg_factors = avg_factors /seasons
    estimation = avg_factors * avg
    estimation_final = estimation[0:n]
    return estimation_final

def expotential_smoothing(int_data,target_data, alpha):
    int_estimates = []
    for i in range(len(int_data)+1):
        if i==0:
            estimate = int_data[0]
        else:
            estimate = alpha * int_data[i-1] + (1-alpha)*int_estimates[-1]
        int_estimates.append(estimate)
    target_estimates = []
    for i in range(len(target_data)+1):
        if i == 0:
            estimate = alpha * int_data[-1] + (1-alpha)*int_estimates[-1]
        else:
            estimate = alpha * target_data[i-1] + (1-alpha)*target_estimates[-1]
        target_estimates.append(estimate)
    return target_estimates
    
def winters_method(int_data,target_data, season_quantity, season_length, forecasting_period):
    forecast_final =np.zeros(forecasting_period)
    avg =np.mean(int_data)
    seasonal_factors = int_data / avg
    seasons = {}
    initial_slope = 0
    avg_factors = np.zeros(season_length)
    for i in range (len(int_data)//season_length):
        avg_factors += seasonal_factors[season_length * i: season_length*i+season_length] 
    avg_factors = avg_factors /season_quantity
    for i in range (season_quantity):
        required_season = int_data[season_length*i : season_length*(i+1) ]
        season_total = np.mean(required_season)
        seasons[i] =season_total
        if i == season_quantity-2:
            initial_slope -= season_total
        elif i == season_quantity-1:
            initial_slope+= season_total
    initial_slope = initial_slope/(season_length**2)
    initial_deseason = seasons[season_quantity-1] + initial_slope * ((season_length-1)/2)
    slope  =  initial_slope
    deseason = initial_deseason
    for i in range (forecasting_period):
        forecast = (deseason + (i+1)*slope)*avg_factors[i%season_length]
        forecast_final[i] += forecast
    return forecast_final

def winters_OSA_parameters(int_data, target_data, season_quantity, season_length, forecasting_period):
    forecast_final = np.zeros(forecasting_period)
    avg =np.mean(int_data)
    seasonal_factors = int_data / avg
    seasons = {}
    initial_slope = 0
    avg_factors = np.zeros(season_length)
    new_factors = avg_factors
    for i in range (len(int_data)//season_length):
        avg_factors += seasonal_factors[season_length * i: season_length*i+season_length] 
    avg_factors = avg_factors /season_quantity
    for i in range (season_quantity):
        required_season = int_data[season_length*i : season_length*(i+1) ]
        season_total = np.mean(required_season)
        seasons[i] =season_total
        if i == season_quantity-2:
            initial_slope -= season_total
        elif i == season_quantity-1:
            initial_slope+= season_total
    initial_slope = initial_slope/(season_length**2)
    initial_deseason = seasons[season_quantity-1] + initial_slope * ((season_length-1)/2)
    slope = initial_slope
    deseason  =  initial_deseason 
    mad_t = 10000000
    for alpha in np.linspace(0,1,11):
       
        for beta in np.linspace(0,1,11):
            
            for gamma in np.linspace(0,1,11):
                
                forecast_final = np.zeros(forecasting_period)
                
                if alpha>= beta:
                    
                    for v in range (forecasting_period):
                        factor = new_factors[v%season_length]
                        demand = target_data[v-1]
                        old_slope = slope
                        old_deseason = deseason
                        deseason = alpha * demand/factor
                        deseason += (1-alpha)* (old_deseason + old_slope)
                        slope = beta * (deseason - old_deseason) + (1-beta)* old_slope
                        new_factors[v%season_length] = gamma * demand / deseason + (1-gamma)*factor
                        forecast = (deseason + slope)* factor
                        forecast_final[v] =forecast      
                    if mad_t >= mad(forecast_final, target_data):
                        mad_t = mad(forecast_final, target_data)
                        best_alpha = alpha
                        best_beta = beta
                        best_gamma = gamma
                    if alpha < 0.3:
                        
                        print(mad(forecast_final, target_data))
                    new_factors =np.zeros(len(new_factors))
                    new_factors += avg_factors
    return (best_alpha, best_beta, best_gamma,)

def winters_OSA(int_data,season_quantity, season_length, forecasting_period, parameters):
    forecast_final = np.zeros(forecasting_period)
    alpha= parameters[0]
    beta= parameters[1]
    gamma= parameters[2]
    avg =np.mean(int_data)
    seasonal_factors = int_data / avg
    seasons = {}
    initial_slope = 0
    avg_factors = np.zeros(season_length)
    for i in range (len(int_data)//season_length):
        avg_factors += seasonal_factors[season_length * i: season_length*i+season_length] 
    avg_factors = avg_factors /season_quantity
    for i in range (season_quantity):
        required_season = int_data[season_length*i : season_length*(i+1) ]
        season_total = np.mean(required_season)
        seasons[i] =season_total
        if i == season_quantity-2:
            initial_slope -= season_total
        elif i == season_quantity-1:
            initial_slope+= season_total
    initial_slope = initial_slope/(season_length**2)
    initial_deseason = seasons[season_quantity-1] + initial_slope * ((season_length-1)/2)
    slope = initial_slope
    deseason  =  initial_deseason 
    for i in range (forecasting_period):
        if i == 0:
            demand = int_data[season_length-1]
        else:  
            demand = forecast_final[i-1]
        factor = avg_factors[i%season_length]
        old_slope = slope
        old_deseason = deseason
        deseason = alpha * demand/factor
        deseason += (1-alpha)* (old_deseason + old_slope)
        slope = beta * (deseason - old_deseason) + (1-beta)* old_slope
        avg_factors[i%season_length] = gamma * demand / deseason + (1-gamma)*factor
        forecast = (deseason + slope)* factor
        forecast_final[i] =forecast

    return forecast_final

def Feasbility(lmd,P):
    ratio = 0
    for i in range (len(lmd)):
        ratio += lmd[i]/P[i]
    print("ratio: " + str(ratio))
    
def T_opt(K, h, lmd, P):
    Ktotal = np.sum(K)
    denominator = 0
    hz= np.zeros(np.size(h))
    for i in range (len(lmd)):
        hz[i] += h[i]*(1-lmd[i]/P[i])
        denominator  += hz[i]*lmd[i]
        print(hz[i]*lmd[i])
    T_opt = math.sqrt(2*Ktotal/(denominator)) 
    return T_opt 

def  T_min(s, P, lmd):
    stotal = np.sum(s)
    Ptotal = np.sum(P)
    lmdtotal = np.sum(lmd)
    T_min = stotal/(1-lmdtotal/Ptotal)
    return T_min

def cost(K, h, lmd, T, P):
    Ktotal = np.sum(K)
    denominator = 0
    hz= np.zeros(np.size(h))
    for i in range (len(lmd)):
        hz[i] += h[i]*(1-lmd[i]/P[i])
        denominator  += hz[i]*lmd[i]
    annual_cost = Ktotal / T +denominator*T/2
    return annual_cost

def time_pct(P,T,s, lmd):
    cycle_prod_time = np.zeros(np.size(P))
    for i in range (len(lmd)):
        cycle_prod_time[i] += lmd[i]*T/P[i]
    total_prod_time = np.sum(cycle_prod_time)/T
    cycle_setup_time= np.sum(s)
    total_setup_time = cycle_setup_time/T
    total_idle_time = 1- total_prod_time- total_setup_time
    prod_pct =   cycle_prod_time* 100/T
    setup_pct = total_setup_time * 100
    idle_pct = total_idle_time *100
    print("prod pct: %" + str(prod_pct))
    print("setup pct:% " + str(setup_pct))
    print("idle pct: %" + str(idle_pct))
    
def QR_Model(lmd, K, p, h, std_dev):
    Q_int = math.sqrt(2*lmd*K/h)
    Rfind_int = Q_int*h/(p*lmd)
    print("Q: " + str(Q_int))
    print("1-F(R) = " + str(Rfind_int))
    print(1-Rfind_int)
    z_int = float(input("Enter z here: "))
    R = z_int*std_dev +lmd/52
    print(R)
    check = input("Do you want to stop?: ")
    while check != "y":
        Lz = float(input("L(z): "))
        nR = std_dev*Lz
        Q = math.sqrt(2*lmd*(K+p*nR)/(h))
        Rfind = Q*h/(p*lmd)
        print("Q: " + str(Q))
        
        print("1-F(R) = " + str(Rfind))
        print(1-Rfind)
        z = float(input("Enter z here: "))
        R = z*std_dev +lmd/52
        print(R)
        check = input("Do you want to stop?: ")
    return  (Q,R)

def QR_Cost(lmd, K, p, h,tau, Variables, std_dev):
    Lz = float(input("Lz: "))
    nR = Lz*std_dev
    cost  = h* (Variables[0]/2+Variables[1] - lmd* tau)+ lmd*K / Variables[0] + p*lmd*nR/Variables[0]
    return cost