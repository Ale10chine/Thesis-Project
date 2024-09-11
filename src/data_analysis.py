import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import time

# Funzione per processare ogni file
def process_file(file_path, fig_num):
    try:
        # Verifica se il file esiste
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Non ho trovato questo file: {file_path}")

        # Carica i dati da un CSV
        data = pd.read_csv(file_path)

        # Pulizia dei dati
        data.replace("/", pd.NA, inplace=True)
        data["PrimalGap"] = pd.to_numeric(data["PrimalGap"], errors='coerce')

        # Filtra i dati per stato "Feasible"
        feasible_data = data[data["Status"] == "Feasible"]

        # Calcola la percentuale di problemi "Feasible" sul totale
        total_problems = len(data)
        feasible_problems = len(feasible_data)
        feasible_percentage = (feasible_problems / total_problems) * 100

        # Calcola il tempo medio per i problemi "Feasible"
        average_time = feasible_data["Time (sec)"].mean()

        # Crea una mappa di colori basata sui valori di "PrimalGap"
        norm = plt.Normalize(feasible_data["PrimalGap"].min(), feasible_data["PrimalGap"].max())
        sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
        sm.set_array([])

        # Visualizzazione del "Optimality Gap" dei problemi "Feasible"
        plt.figure(figsize=(16, 9), num=f'Figura {fig_num}')  # Imposta una dimensione grande (es. 16:9)
        ax = sns.barplot(x='ProblemName', y='PrimalGap', data=feasible_data, palette="YlOrRd", hue='PrimalGap', dodge=False, legend=False)
        plt.xlabel('')

        plt.title(f'Primal Gap dei Problemi Feasible - {file_path}')
        plt.ylabel('Primal Gap')

        # Migliora l'aspetto delle etichette dell'asse X
        plt.xticks(rotation=45, ha='right')

        # Aggiungi la barra dei colori
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
        cbar.set_label('Primal Gap')

        # Aggiungi il tempo medio e le considerazioni sotto il grafico
        info_text = (
            f'Percentuale di Problemi Feasible: {feasible_percentage:.2f}%\n'
            f'Tempo Medio per Trovare una Soluzione: {average_time:.2f} sec\n\n'
        )

        # Regola lo spazio sotto il grafico
        plt.subplots_adjust(bottom=0.30)

        # Aggiungi il testo sotto il grafico
        plt.figtext(0.1, 0.1, info_text, ha='left', va='top', fontsize=10, color='black', bbox=dict(facecolor='lightgrey', alpha=0.5))

        # Salva il grafico nella cartella "../charts" con un nome file basato sul numero della figura
        chart_path = f"../charts/Figura_{fig_num}.png"
        os.makedirs(os.path.dirname(chart_path), exist_ok=True)  # Crea la cartella se non esiste
        plt.savefig(chart_path, dpi=300, bbox_inches='tight')  # Salva a alta risoluzione con bounding box "tight"

        
        plt.show()
        # Chiudi la finestra del grafico
        plt.close()

    except FileNotFoundError as e:
        print(e)

# Lista dei file da processare
file_paths = [
    '../out_csv/seed1_p1/Result.csv',
    '../out_csv/seed1_p2/Result.csv',
    '../out_csv/seed1_p3/Result.csv'
]

# Esegui la funzione per ogni file
for i, file_path in enumerate(file_paths, start=1):
    process_file(file_path, fig_num=i)
