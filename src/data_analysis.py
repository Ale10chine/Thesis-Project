import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import time

# Function to iter in every file
def process_file(file_path, fig_num):
    try:
        # Check if file exists 
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Non ho trovato questo file: {file_path}")

        # Load data from a csv file 
        data = pd.read_csv(file_path)

        # Cleaning of data
        data.replace("/", pd.NA, inplace=True)
        data["PrimalGap"] = pd.to_numeric(data["PrimalGap"], errors='coerce')

        # Look for problems that are "Feasible"
        feasible_data = data[data["Status"] == "Feasible"]

        # Calculate the percentage of "Feasible" problems on the total
        total_problems = len(data)
        feasible_problems = len(feasible_data)
        feasible_percentage = (feasible_problems / total_problems) * 100

        # Calculate the average time to find solution for "Feasible" problems
        average_time = feasible_data["Time (sec)"].mean()

        # Calculate the average Primal Gap for "Feasible" problems 
        average_primal_gap = feasible_data["PrimalGap"].mean()

        # Sorts the data by ascending primal gap
        feasible_data = feasible_data.sort_values(by="PrimalGap", ascending=True)

        # Create a color map based on the values ​​of "PrimalGap"
        norm = plt.Normalize(feasible_data["PrimalGap"].min(), feasible_data["PrimalGap"].max())
        sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
        sm.set_array([])

        # Plot the charts
        plt.figure(figsize=(40, 15), num=f'Figura {fig_num}')  # For size of chart
        ax = sns.barplot(x='ProblemName', y='PrimalGap', data=feasible_data, palette="YlOrRd", hue='PrimalGap', dodge=False, legend=False)
        plt.xlabel('')

        plt.title(f'Primal Gap dei Problemi Feasible - {file_path}')
        plt.ylabel('Primal Gap')

        # Rotate the label of problems names
        plt.xticks(rotation=45, ha='right')

        # Add the color bar
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
        cbar.set_label('Primal Gap')

        # Add stats under the chart
        info_text = (
            f'Problemi Feasible: {feasible_problems} su {total_problems} ({feasible_percentage:.2f}%)\n'
            f'Tempo Medio per Trovare una Soluzione: {average_time:.2f} sec\n'
            f'Media del Primal Gap (Feasible): {average_primal_gap:.4f}\n\n'
        )

        # Adjuste space under the chart
        plt.subplots_adjust(bottom=0.30)

        # Adjuste the stats box under the chart
        plt.figtext(0.1, 0.1, info_text, ha='left', va='top', fontsize=10, color='black', bbox=dict(facecolor='lightgrey', alpha=0.5))

        # Save the chart to the "../charts" folder with a file name based on the figure number
        chart_path = f"../charts/Figura_{fig_num}.png"
        os.makedirs(os.path.dirname(chart_path), exist_ok=True)  # Create directory if doesn't exists
        plt.savefig(chart_path, dpi=300, bbox_inches='tight')  # Save high resolution with "tight" bounding box

        plt.show()
        # Close the window of the chart
        plt.close()

    except FileNotFoundError as e:
        print(e)

# List of file to execute
file_paths = [


    '../out_csv/seed1_p1/1.csv',
    '../out_csv/seed1_p1/I5_T30_s1.csv',
    '../out_csv/seed1_p1/I15_T30_s1.csv',
    '../out_csv/seed1_p1/I20_T30_s1.csv',
    '../out_csv/seed1_p1/T15_s1.csv',
    '../out_csv/seed1_p1/T25_s1.csv',
    '../out_csv/seed1_p1/T35_s1.csv',
    '../out_csv/seed1_p1/T45_s1.csv',

    '../out_csv/seed1_p2/2.csv',
    
    '../out_csv/seed1_p3/3.csv',

    '../out_csv/seed2_p1/4.csv',
    '../out_csv/seed2_p1/I5_T30_s2.csv',
    '../out_csv/seed2_p1/I15_T30_s2.csv',
    '../out_csv/seed2_p1/I20_T30_s2.csv',
    '../out_csv/seed2_p1/T15_s2.csv',
    '../out_csv/seed2_p1/T25_s2.csv',
    '../out_csv/seed2_p1/T35_s2.csv',
    '../out_csv/seed2_p1/T45_s2.csv',

    '../out_csv/seed2_p2/5.csv',
    
    '../out_csv/seed2_p3/6.csv',

    '../out_csv/seed3_p1/7.csv',
    '../out_csv/seed3_p1/I5_T30_s3.csv',
    '../out_csv/seed3_p1/I15_T30_s3.csv',
    '../out_csv/seed3_p1/I20_T30_s3.csv',
    '../out_csv/seed3_p1/T15_s3.csv',
    '../out_csv/seed3_p1/T25_s3.csv',
    '../out_csv/seed3_p1/T35_s3.csv',
    '../out_csv/seed3_p1/T45_s3.csv',

    '../out_csv/seed3_p2/8.csv',
    
    '../out_csv/seed3_p3/9.csv',


 

]
    #'../out_csv/seed1_p1/I20_T30_s1.csv',
    #'../out_csv/seed2_p1/I20_T30_s2.csv',
    #'../out_csv/seed3_p1/I20_T30_s3.csv',

    #'../T15_s1.csv',
    #'../T15_s2.csv',
    #'../T15_s3.csv',



    #'../out_csv/seed1_p1/1.csv',
    #'../2.csv',
    #'../3.csv',
    #'../T25_s1.csv',
    #'../T25_s2.csv',
    #'../T25_s3.csv',
    #'../T35_s1.csv',
    #'../T35_s2.csv',
    #'../T35_s3.csv',
    #'../T45_s1.csv',
    #'../T45_s2.csv',
    #'../T45_s3.csv',
    #'../I15_T30_s1.csv',
    #'../I15_T30_s2.csv',
    #'../I15_T30_s3.csv'



#    '../T25_s1.csv',
#    '../T25_s2.csv',
#    '../T25_s3.csv',
#    '../T35_s1.csv',
#    '../T35_s2.csv',
#    '../T35_s3.csv',
#    '../T45_s1.csv',
#    '../T45_s2.csv',
#    '../T45_s3.csv'
#    '../I15_T30_s1.csv',
#    '../I15_T30_s2.csv',
#    '../I15_T30_s3.csv'





#    '../out_csv/seed1_p1/1.csv',
#    '../out_csv/seed1_p2/2.csv',
#    '../out_csv/seed1_p3/3.csv',
#    '../out_csv/seed2_p1/4.csv',
#    '../out_csv/seed2_p2/5.csv',
#    '../out_csv/seed2_p3/6.csv',
#    '../out_csv/seed3_p1/7.csv',
#    '../out_csv/seed3_p2/8.csv',
#    '../out_csv/seed3_p3/9.csv',

#'../T35_s3.csv',

# Iteration of the function for every file
for i, file_path in enumerate(file_paths, start=1):
    process_file(file_path, fig_num=i)
