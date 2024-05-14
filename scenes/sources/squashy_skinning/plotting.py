import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output_kappa.csv')
kappa = np.array(df['kappa'])
time = np.array([i for i, value in enumerate(df['kappa'])])

# Defined the range of k values for which we want to plot the graph
# k_values = np.arange(0, 50)

# Calculate corresponding kappa values using the scene_model::kappa function
# kappa_values = [scene_model.kappa(k) for k in k_values]

plt.plot(time, kappa)
plt.xlabel( 'time t    (0.017*k)' )
plt.ylabel( 'kappa_function' )
plt.title( 'Graph of scene_model::kappa' )
plt.grid( True )
plt.show()