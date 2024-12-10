"""
Nombre: Circuito RLC en Serie
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

def F3(ecuacion,t,i,u,R,L,C):
     g = 9.81
     pi = 3.141592653589793238  
     exp = 2.71828182845904523
     return eval(str(ecuacion))

def Grafica(Ft,tf):
     h = np.linspace(-100,100,2)
     g = []
     for i in h:
          g.append(0)

     archivo = open("Circuito RLC Serie Metodo Euler.txt","w")
     archivo.write("t , It\n")
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
     ax.set_ylabel("Corriente (A)")
     #plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0))
     #plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
     ax.plot(h,g,color="red")
     ax.plot(g,h,color="red")
     ax.plot(0,I,color="white",marker =".")
     ax.plot(tf,Ft,color="yellow")#,label=f"$f(x) = {Vt}$")
     print(type(len(tf)-1))
     ax.set_xlim(0,tf[len(tf)-1])
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

R = 10 # Ohms
L = 1e-3 # Henrios
C = 35e-6 # Faradios
I = 0 # Amperios -- Corriente inicial del inductor
V = -75e-3 # Voltios -- Voltaje inicial del capacitor
tf = 1.5e-3 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2
h = tf/1000

vL = -((I*R)+V) # Voltaje en el inductor

t = 0
i = I
u = vL/L
f = " -(R/L)*u - (1/(L*C))*i "
il = []
ul = []
tl = []

while t < tf:
     
     if t == 0:
          il.append(i)
          ul.append(u)
          tl.append(t)
     
     i = i + h*u
     u = u + h*(F3(f,t,i,u,R,L,C))
     il.append(i)
     ul.append(u)
     tl.append(t)
     
     t += h

Grafica(il,tl)
