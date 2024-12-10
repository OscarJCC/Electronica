"""
Nombre: Circuito RLC en Paralelo
Fecha Inicio:
Fecha Final:
Autor: Oscar Joel Castro Contreras y David Peralez
Descripcion:

"""

def F(ecuacion,y):
     g = 9.81
     pi = 3.141592653589793238  
     exp = 2.71828182845904523
     return eval(str(ecuacion))

def F3(ecuacion,t,v,u,R,L,C):
     g = 9.81
     pi = 3.141592653589793238  
     exp = 2.71828182845904523
     return eval(str(ecuacion))

def Grafica(Ft,tf):
     k = np.linspace(-100,100,2)
     g = []
     for i in k:
          g.append(0)

     archivo = open("Circuito RLC Paralelo Metodo Euler.txt","w")
     archivo.write("t , Vt\n")
     A = len(Ft)
     for i in range(A):
          archivo.write(f"{tf[i]} , ")
          archivo.write(f"{Ft[i]}\n")
     archivo.close()
     
     plt.style.use("dark_background") 
     fig, ax = plt.subplots()
     ax.grid(color="green")
     #ax.set_title(N)
     ax.set_xlabel("Tiempo (s)")
     ax.set_ylabel("Voltaje (V)")
     #plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0))
     #plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
     ax.plot(k,g,color="red")
     ax.plot(g,k,color="red")
     ax.plot(0,V,color="white",marker =".")
     ax.plot(tf,Ft,color="yellow")#,label=f"$f(x) = {Vt}$")
     ax.set_xlim(-tf[len(tf)-1]/100,tf[len(tf)-1])
     limy = float(OrCre(Ft))
     ax.set_ylim(-limy,limy)
     plt.show()

def OrCre(lista):
     for i in range(1,len(lista)):
          for j in range(len(lista)-i):
               if(lista[j] > lista[j+1]):
                    aux = lista[j]
                    lista[j] = lista[j+1]
                    lista[j+1] = aux
     
     if abs(lista[0]) >= lista[len(lista)-1]:
          return abs(lista[0])
     else:
          return abs(lista[len(lista)-1])

import numpy as np
import matplotlib.pyplot as plt

R = 150 # Ohms
L = 70e-3 # Henrios
C = 70e-6 # Faradios
I = 200e-3 # Amperios -- Corriente inicial del inductor
V = -15 # Voltios -- Voltaje inicial del capacitor
tf = .1 # Tiempo final(segundos)
h = tf/1000

ic = -((abs(V/R))+I)

t = 0
v = V
u = ic/C
f = " - (1/(R*C))*u - (1/(L*C))*v "
vl = []
ul = []
tl = []
while t < tf:
     
     if t == 0:
          vl.append(v)
          ul.append(u)
          tl.append(t)

     v = v + h*u
     u = u + h*(F3(f,t,v,u,R,L,C))

     vl.append(v)
     ul.append(u)
     tl.append(t)
     
     t += h

Grafica(vl,tl)
