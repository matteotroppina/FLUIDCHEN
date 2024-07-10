import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

csv_files = {
    "0.0001": "iterations_nu_0.0001.csv",
    "0.0005": "iterations_nu_0.0005.csv",
    "0.002": "iterations_nu_0.002.csv",
    "0.01": "iterations_nu_0.01.csv"
}

total_time = 50
plt.figure()

for nu, filename in csv_files.items():
    if not os.path.exists(filename):
        print(f"File {filename} does not exist. Skipping.")
        continue

    df = pd.read_csv(filename, header=None, sep=",").T

    n_steps = len(df)
    timestep = np.linspace(0, total_time, n_steps)

    plt.scatter(timestep, df, label=f'nu={nu}', s=1)

    plt.title(f'Iterations vs. Time for different nu, SOR omg=1.7')
    plt.xlabel('Time')
    plt.ylabel('Number of Iterations')

plt.grid()
plt.legend(loc="lower right")
plt.savefig(f'nu_iterations.png')

