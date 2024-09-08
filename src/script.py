import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Carica i dati da un CSV
data = pd.read_csv('../out_csv/Result.csv')

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
plt.figure(figsize=(8, 4))  # Riduci la dimensione della figura per evitare che sia troppo grande
ax = sns.barplot(x='ProblemName', y='PrimalGap', data=feasible_data, palette="YlOrRd", hue='PrimalGap', dodge=False, legend=False)
plt.xlabel('')

plt.title('Primal Gap dei Problemi in cui è stata trovata una sluzione ammissibile con ACS')
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
plt.subplots_adjust(bottom=0.30)  # Aumenta lo spazio sotto il grafico

# Aggiungi il testo sotto il grafico
plt.figtext(0.1, 0.1, info_text, ha='left', va='top', fontsize=10, color='black', bbox=dict(facecolor='lightgrey', alpha=0.5))

plt.show()
plt.close()
