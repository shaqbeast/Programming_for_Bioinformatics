#!/usr/bin/env python3
import sys 
import csv 
import matplotlib.pyplot as plt 

covid_file = sys.argv[1]
year = [] 
anomaly = []

with open(covid_file, 'r') as file: 
    # read the covid file as a csv file
    reader_csv = csv.reader(file)
    
    # skip the first 5  lines of the csv file
    for i in range(5):
        next(reader_csv)
    
    # save the values of the csv columns into year and anomaly lists respectively 
    for row in reader_csv:
        year.append(int(row[0]))
        anomaly.append(float(row[1]))

# plot the year as x-axis and anomaly values as the y-axis 
plt.plot(year, anomaly) 

# add titles and labels
plt.title('Global Average Temperature Anomalies')
plt.xlabel('Year')
plt.ylabel('Anomaly')

# add a grid for better visualization
plt.grid()
plt.show()
        
    