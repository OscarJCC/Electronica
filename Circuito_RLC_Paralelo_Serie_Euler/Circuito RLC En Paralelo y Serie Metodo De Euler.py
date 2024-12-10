"""
Nombre: Circuito RLC En Paralelo y Serie Metodo De Euler
Fecha Inicio:
Fecha Final:
Autor: Oscar Joel Castro Contreras y David Peralez
Descripcion:

"""
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def Grafica(N,S,Ft,tf):
     k = np.linspace(-100,100,2)
     g = []
     for i in k:
          g.append(0)
     
     print(g)

     archivo = open("Circuito RLC En Paralelo y Serie Metodo de Euler.txt","w")
     if S == 1:
          archivo.write("t , Vt\n")
     elif S == 2:
          archivo.write("t , It\n")
     A = len(Ft)
     for i in range(A):
          archivo.write(f"{tf[i]} , ")
          archivo.write(f"{Ft[i]}\n")
     archivo.close()
     
     plt.style.use("dark_background") 
     fig, ax = plt.subplots()
     ax.grid(color="green")
     ax.set_title(N)
     ax.set_xlabel("Tiempo (s)")
     #plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0))
     #plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
     ax.plot(k,g,color="red")
     ax.plot(g,k,color="red")
     ax.plot(tf,Ft,color="yellow")#,label=f"$f(x) = {Vt}$")
     if S == 1:
          ax.plot(0,V,color="white",marker =".")
          ax.set_ylabel("Voltaje (V)")
     elif S == 2:
          ax.plot(0,I,color="white",marker =".")
          ax.set_ylabel("Corriente (A)")
     ax.set_xlim(-tf[len(tf)-1]/100,tf[len(tf)-1])
     limy = float(OrCre(Ft))
     ax.set_ylim(-limy,limy)
     plt.savefig("Circuito RLC En Paralelo y Serie Metodo de Euler.pdf")
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

def F3(ecuacion,t,l,u,R,L,C):
     return eval(str(ecuacion))

def Euler(S,R,L,C,l,f,u,h):
     t = 0
     lf = []
     tl = []
     while t < tf:
          if t == 0:
               lf.append(l)
               tl.append(t)
          l = l + h*u
          u = u + h*(F3(f,t,l,u,R,L,C))
          lf.append(l)
          t += h
          tl.append(t)
     Grafica(N,S,lf,tl)
     
R = 150 # Ohms
L = 70e-3 # Henrios
C = 70e-6 # Faradios
I = 200e-3 # Amperios -- Corriente inicial del inductor
V = -15 # Voltios -- Voltaje inicial del capacitor
tf = .1 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2
h = tf/1000

if S == 1: # Paralelo

     a = 1/(2*R*C) # Alpha
     w = 1/(np.sqrt(L*C)) # w_o
     
     ic = -((-(V/R))+I) # Corriente en el capacitor
     
     l = v = V
     
     u = ic/C
     
     f = " - (1/(R*C))*u - (1/(L*C))*l "
     
     if a > w:

          N = "Circuito RLC En Paralelo Metodo de Euler: Sistema Sobre Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,v,f,u,h)

     elif a == w:

          N = "Circuito RLC En Paralelo Metodo de Euler: Sistema Criticamente Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,v,f,u,h)

     elif a < w:

          N = "Circuito RLC En Paralelo Metodo de Euler: Sistema Sub Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,v,f,u,h)

elif S == 2: #Serie
     
     a = R/(2*L) # Alpha
     w = 1/np.sqrt(L*C) # w_o
     
     vL = -((I*R)+V) # Voltaje en el inductor
     
     l = i = I
     
     u = vL/L
     
     f = " -(R/L)*u - (1/(L*C))*l "
     
     if a > w:
          
          N = "Circuito RLC En Serie Metodo de Euler: Sistema Sobre Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,i,f,u,h)
          
     elif a == w:
          
          N = "Circuito RLC En Serie Metodo de Euler: Sistema Criticamente Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,i,f,u,h)
               
     elif a < w:
          
          N = "Circuito RLC En Serie Metodo de Euler: Sistema Sub Amortiguado"
          print("\n\n\t",N)
          
          Euler(S,R,L,C,i,f,u,h)
     
"""
-----Circuito RLC En Paralelo: -----

     Sistema Sobre Amortiguado:
     
R = 20 # Ohms
L = 25e-3 # Henrios
C = 10e-6 # Faradios
I = 10e-3 # Amperios -- Corriente inicial del inductor
V = 0 # Voltios -- Voltaje inicial del capacitor
tf = .01 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2
          
          
     Sistema Criticamente Amortiguado:
     
R = 20 # Ohms
L = 16e-3 # Henrios
C = 10e-6 # Faradios
I = 3 # Amperios -- Corriente inicial del inductor
V = -2.5 # Voltios -- Voltaje inicial del capacitor
tf = .004 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2


     Sistema Sub Amortiguado:
     
R = 150 # Ohms
L = 70e-3 # Henrios
C = 70e-6 # Faradios
I = 200e-3 # Amperios -- Corriente inicial del inductor
V = -15 # Voltios -- Voltaje inicial del capacitor
tf = .1 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2


-----Circuito RLC En Serie: -----

     Sistema Sobre Amortiguado:
     
R = 1e3 # Ohms
L = 3e-3 # Henrios
C = 1.2e-3 # Faradios
I = 300e-3 # Amperios -- Corriente inicial del inductor
V = 5 # Voltios -- Voltaje inicial del capacitor
tf = .00003 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2
          
          
     Sistema Criticamente Amortiguado:
     
R = 100 # Ohms
L = 25e-3 # Henrios
C = 10e-6 # Faradios
I = 10 # Amperios -- Corriente inicial del inductor
V = 0 # Voltios -- Voltaje inicial del capacitor
tf = .004 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2


     Sistema Sub Amortiguado:
     
R = 10 # Ohms
L = 70e-3 # Henrios
C = 70e-6 # Faradios
I = 200e-3 # Amperios -- Corriente inicial del inductor
V = -15 # Voltios -- Voltaje inicial del capacitor
tf = .1 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2
"""