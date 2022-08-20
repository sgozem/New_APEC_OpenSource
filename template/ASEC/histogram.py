import numpy as np
import matplotlib.pyplot as plt

name = 'NAME'  # put the name label here
f = open('Energies', 'r')
a = f.read()
f.close()

a = a.split()
a = a[10::4]
a = [float(i) for i in a]

#plt.hist(a)

#options: http://matplotlib.org/api/pyplot_api.html

#bins = np.linspace(50,58,20)i  # init,final,num_intervals
#plt.hist(a, bins=20)
plt.hist(a, bins=15, alpha=0.8, color='red', edgecolor='black', width=0.20)
plt.title(name)
plt.xlabel('Abs Energy')
plt.ylabel('Freq.')
plt.savefig('%s.png' % name)
#plt.show()



