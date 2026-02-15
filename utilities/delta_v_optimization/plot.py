import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

data = pd.read_csv(Path(__file__).parent / 'delta_v_optimization.csv')

# Filter out inf values so the plot is readable
data = data[np.isfinite(data['delta_v'])]

orbit_samples = data['orbit_samples']
delta_v = data['delta_v']

plt.scatter(orbit_samples, delta_v)
plt.xlabel('Orbit Samples')
plt.ylabel('Delta V (m/s)')
plt.title('Delta V Optimization')
plt.show()