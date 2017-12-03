#este script lee los datos producidos por el programa en C y los grafica.
import numpy as np
import matplotlib.pyplot as plt

x= np.linspace(0,300,300)
Datos = np.genfromtxt('new_data.dat',delimiter= ' ')
Y =Datos[:,0]
Y2=Datos[:,1]




plt.plot(x,Y,label="Velocidad observada")
plt.plot(x,Y2,label="Velocidad del modelo")
plt.xlabel('enumeracion de los datos')
plt.ylabel('datos')
plt.title('grafica de comparacion')
plt.legend()
plt.tight_layout()


plt.savefig('grafica.png')
#plt.show()
plt.close()


