#!/usr/bin/env python

''' 
This script plots and calculates circular-linear correlations between
all possible combinations of circular variables x and linear variables y 

By Luis Sosa
'''

from mpmath import *
import random
import sys
from pandas import Series
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from scipy import mean
from statistics import median
from copy import copy
import csv

def calculate_correlation_coef(cicrcular, linear, parametric=True):
    ''' Calculates correlation coefficient '''

    # Get dataset length
    assert (len(linear) == len(circular)), "The length of the circular and linear variable arrays must match"
    n = len(linear)

    # Parametric test
    if parametric:
        c = [float(cos(entry)) for entry in circular]
        s = [float(sin(entry)) for entry in circular]
        rxc = pearsonr(linear, c)[0]
        rxs = pearsonr(linear, s)[0]
        rcs = pearsonr(c, s)[0]
        r2 = (rxc**2 + rxs**2 - 2*rxc*rxs*rcs)/(1-rcs**2)
        return sqrt(r2)
    
    # Non-parametric test
    else:
        # Ranked data
        cir_ranks = Series(circular).rank()
        lin_ranks = Series(linear).rank()
        
        # Interim calculations
        beta = [2*pi*entry/n for entry in cir_ranks]
        C = sum([x_i*cos(beta_i) for beta_i,x_i in zip(beta, lin_ranks)])
        S = sum([x_i*sin(beta_i) for beta_i,x_i in zip(beta, lin_ranks)])

        # Denormalized correlation coefficient
        # U = 24*(C**2 + S**2)/(n**2 * (n+1))

        # Normalization coefficient
        alpha = 2*(sin(pi/n))**4 / (1 + (cos(pi/n)))**3 if n%2 \
                else 1/(1 + 5*(cot(pi/n))**2 + 4*(cot(pi/n))**4) 

        # Normalized correlation coefficient
        D = alpha*(C**2 + S**2)
        return D


def calculate_correlation(circular, linear, parametric=True, iterations=100):
    ''' Calculates circular statistics and performs permutation test to determine p-value '''
    actual_correl = calculate_correlation_coef(circular, linear, parametric=parametric)
    stronger_correls = 0
    rand_linear = copy(linear)

    for i in range(iterations):
        # Continually shuffle data and recalculate correlation
        random.shuffle(rand_linear)
        simulated_correl = calculate_correlation_coef(circular, rand_linear, parametric=parametric)

        # Note random correlation as strong or stronger than observed correlation
        if abs(simulated_correl) >= actual_correl:
            stronger_correls += 1
    
    # p-value is the ratio of stronger random correlations there are
    p_val = stronger_correls/iterations
    return (round(actual_correl,3), round(p_val,3))

def plot_data(circular, linear, correl, pval, title="", dual=True):  
    ''' Plots data and correlation statistics on circular histogram '''

    # Create plot
    plt.figure()
    ax = plt.subplot(111, polar=True)
    bars = ax.bar(circular, linear, width=2*pi/len(circular), edgecolor='black')
    
    # Use custom colors and opacity
    for r, bar in zip(linear, bars):
        bar.set_facecolor((random.random(), 1, random.random()))
        bar.set_alpha(0.8)
    
    # Adjust labelling
    plt.title(f"{title}\n")
    ax.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])
    ax.set_theta_offset(pi/2)
    ax.set_theta_direction(-1)

    # Add corelation data
    if dual:
        plt.figtext(0.15, 0.1, f"Parametric:\nR = {correl[0]}\np-val={p_val[0]}")
        plt.text(0.9, 0, f"Non-parametric:\nR = {correl[1]}\np-val={p_val[1]}", transform=ax.transAxes)
    else:
        plt.figtext(0.15, 0.1, f"R = {correl}\np-val={p_val}")
    
    # Show plot non-blockingly
    plt.ion()
    plt.show()
    
def read_input(datafile, degrees=True):
    with open(datafile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        linear_columns = []
        circular_columns = []
        compounded_series = []
        
        linear_names = []
        circular_names = []
        
        for row in csv_reader:

            # Detect and store series header
            if line_count == 0:
                for i in range(len(row)):
                    # Headers preceeded by _x_ are circular series
                    if row[i].find("_x_") == 0:
                        circular_columns.append(i)
                        circular_names.append(row[i][3:])

                    # Headers preceeded by _y_ are linear series
                    elif row[i].find("_y_") == 0:
                        linear_columns.append(i)
                        linear_names.append(row[i][3:])

                # Create a list of tuples to accomodate all pairings of linear and circular variables
                compounded_series = [[[],[]] for _ in range((len(linear_columns)*len(circular_columns)))]
            
            # Store data with its corresponding series
            else:
                for i in range(len(circular_columns)):
                    circular = circular_columns[i]
                    # Exclude missing data
                    if row[circular] == "-1":
                        continue

                    for j in range(len(linear_columns)):
                        linear = linear_columns[j]
                        # Exclude missing data
                        if row[linear] == "-1":
                            continue

                        # Compounded series will contain all combinations of linear and circular variables
                        else:
                            x = float(row[circular])
                            y = float(row[linear])

                            # Cast angular data into radians
                            if degrees:
                                x = float(x*pi/180)
                                pass

                            compounded_series[i*len(linear_columns) + j][0].append(x)
                            compounded_series[i*len(linear_columns) + j][1].append(y)
                            
            line_count += 1

        print(f'Processed {line_count} lines, producing {len(compounded_series)} distinct combinations')

        titles = []
        for x in circular_names:
            for y in linear_names:
                titles.append(f'{x} vs {y}')

        return compounded_series, titles
        
    

if __name__ == "__main__":
    print('ass')
    quit()
    # Get input file location
    if len(sys.argv) < 2:
        datafile = input("Please input the filepath of the .csv file you wish to analyze: ")
    else:
        datafile = sys.argv[1]

    # Relevance contains the minimum correlation and p-value that are considered relevant enough to plot
    relevance = (0.2, 0.1)

    # Parse input and generate compounded data
    compounded_data, titles = read_input(datafile)

    # Calculate correlation and plot data for each circular-linear pair
    for i in range(len(compounded_data)):
        series = compounded_data[i]
        circular = series[0]
        linear = series[1]
        correl = [None] * 2
        p_val = [None] * 2

        correl[0], p_val[0] = calculate_correlation(circular, linear, parametric=True, iterations=1000)
        correl[1], p_val[1] = calculate_correlation(circular, linear, parametric=False, iterations=1000)
        
        # Display only relevant correlations
        if( (correl[0] >= relevance[0] and p_val[0] <= relevance[1]) or (correl[1] >= relevance[0] and p_val[1] <= relevance[1]) ):
            plot_data(circular, linear, correl, p_val, title=titles[i], dual=True)

        print(f"{round((i+1)/len(compounded_data) * 100,0)}% complete")
       
    input("Press any key to exit...")    