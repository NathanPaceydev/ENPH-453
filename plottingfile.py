#please work, coding gods
import numpy as np 
import matplotlib.pyplot as plt
  
channel= np.loadtxt("pulse100xmul.CNF.txt", usecols=1, dtype=float) 
energy= np.loadtxt("pulse100xmul.CNF.txt", usecols=2, dtype=float) 
counts= np.loadtxt("pulse100xmul.CNF.txt", usecols=3, dtype=float) 
#rate= np.loadtxt("pulse100xmul.CNF.txt", usecols=4, dtype=float) 

plt.plot(channel, energy)
plt.show()