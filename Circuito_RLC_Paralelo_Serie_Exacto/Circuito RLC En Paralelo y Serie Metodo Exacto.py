"""
Nombre: Circuito RLC En Paralelo y Serie Metodo Exacto
Fecha Inicio:
Fecha Final:
Autor: Oscar Joel Castro Contreras y David Peralez
Descripcion:

"""
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

def Grafica(N,S,Ft,tf):
     k = np.linspace(-100,100,2)
     g = []
     for i in k:
          g.append(0)

     c = np.linspace(0,tf,1001)
     d = []
     for i in c:
          e = Ft.subs(t,i)
          d.append(e)

     archivo = open("Circuito RLC En Paralelo y Serie Metodo Exacto.txt","w")
     archivo.write("t , It\n")
     A = len(d)
     for i in range(A):
          archivo.write(f"{c[i]} , ")
          archivo.write(f"{d[i]}\n")
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
     ax.plot(c,d,color="lightblue")#,label=f"$f(x) = {Vt}$")
     if S == 1:
          ax.plot(0,V,color="white",marker =".")
          ax.set_ylabel("Voltaje (V)")
     elif S == 2:
          ax.plot(0,I,color="white",marker =".")
          ax.set_ylabel("Corriente (A)")
     ax.set_xlim(-tf/100,tf)
     limy = float(OrCre(d))
     ax.set_ylim(-limy,limy)
     plt.savefig("Circuito RLC En Paralelo y Serie Metodo Exacto.pdf")
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

R = 10 # Ohms
L = 1e-3 # Henrios
C = 35e-6 # Faradios
I = 0 # Amperios -- Corriente inicial del inductor
V = -75e-3 # Voltios -- Voltaje inicial del capacitor
tf = 1.5e-3 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2

if S == 1: # Paralelo
     
     t,A1,A2 = symbols("t,A1,A2")

     a = 1/(2*R*C) # Alpha
     w = 1/(np.sqrt(L*C)) # w_o
     
     ic = -((-(V/R))+I) # Corriente en el capacitor

     if a > w:

          N = "Circuito RLC En Paralelo Metodo Exacto: Sistema Sobreamortiguado"
          print("\n\n\t",N)

          m1 = -a + sqrt(a**2-w**2)
          m2 = -a - sqrt(a**2-w**2)

          sol = list(linsolve([A1 + A2 - V,m1*A1 + m2*A2 - ic/C],(A1,A2)))
          A1 = sol[0][0]
          A2 = sol[0][1] 
          
          Vt = A1*exp(m1*t) + A2*exp(m2*t)
          print("\nVt = ",f"{Vt}")

     elif a == w:

          N = "Circuito RLC En Paralelo Metodo Exacto: Sistema Criticamente amoritiguado"
          print("\n\n\t",N)

          A1 = V
          A2 = (ic/C)-(-a*V)

          Vt = exp(-a*t)*(A1 + A2*t)
          print("\nVt = ",f"{Vt}")

     elif a < w:

          N = "Circuito RLC En Paralelo Metodo Exacto: Sistema Subamoritiguado"
          print("\n\n\t",N)

          wd = sqrt(w**2-a**2)

          B1 = V
          B2 = ((ic/C)-(-a*V))/wd

          Vt = exp(-a*t)*(B1*cos(wd*t) + B2*sin(wd*t))
          print("\nVt = ",f"{Vt}")
     
     Grafica(N,S,Vt,tf)

elif S == 2: #Serie
     
     t,A1,A2 = symbols("t,A1,A2")
     
     a = R/(2*L) # Alpha
     w = 1/sqrt(L*C) # w_o
     
     vL = -((I*R)+V) # Voltaje en el inductor
     
     if a > w:
          
          N = "Circuito RLC En Serie Metodo Exacto: Sistema Sobreamortiguado"
          print("\n\n\t",N)

          m1 = -a + sqrt(a**2-w**2)
          m2 = -a - sqrt(a**2-w**2)

          sol = list(linsolve([A1 + A2 - I,m1*A1 + m2*A2 - vL/L],(A1,A2)))
          A1 = sol[0][0]
          A2 = sol[0][1] 

          It = A1*exp(m1*t) + A2*exp(m2*t)
          print("\nIt = ",f"{It}")
          
     elif a == w:
          
          N = "Circuito RLC En Serie Metodo Exacto: Sistema Criticamente amoritiguado"
          print("\n\n\t",N)

          A1 = I
          A2 = (vL/L)-(-a*I)

          It = exp(-a*t)*(A1 + A2*t)
          print("\nIt = ",f"{It}")
               
     elif a < w:
          
          N = "Circuito RLC En Serie Metodo Exacto: Sistema Subamoritiguado"
          print("\n\n\t",N)

          wd = sqrt(w**2-a**2)

          B1 = I
          B2 = ((vL/L)-(-a*I))/wd

          It = exp(-a*t)*(B1*cos(wd*t) + B2*sin(wd*t))
          print("\nIt = ",f"{It}")
     
     Grafica(N,S,It,tf)

"""
-----Circuito RLC En Paralelo:-----

     Sistema Sobreamortiguado:
     
R = 20 # Ohms
L = 25e-3 # Henrios
C = 10e-6 # Faradios
I = 10e-3 # Amperios -- Corriente inicial del inductor
V = 0 # Voltios -- Voltaje inicial del capacitor
tf = .01 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2
          
          
     Sistema Criticamente amoritiguado:
     
R = 20 # Ohms
L = 16e-3 # Henrios
C = 10e-6 # Faradios
I = 3 # Amperios -- Corriente inicial del inductor
V = -2.5 # Voltios -- Voltaje inicial del capacitor
tf = .004 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2


     Sistema Subamoritiguado:
     
R = 150 # Ohms
L = 70e-3 # Henrios
C = 70e-6 # Faradios
I = 200e-3 # Amperios -- Corriente inicial del inductor
V = -15 # Voltios -- Voltaje inicial del capacitor
tf = .1 # Tiempo final(segundos)
S = 1 # Circuto en Paralelo = 1, Serie = 2


-----Circuito RLC En Serie-----

     Sistema Sobreamortiguado:
     
R = 1e3 # Ohms
L = 3e-3 # Henrios
C = 1.2e-3 # Faradios
I = 300e-3 # Amperios -- Corriente inicial del inductor
V = 5 # Voltios -- Voltaje inicial del capacitor
tf = .00003 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2
          
          
     Sistema Criticamente amoritiguado:
     
R = 100 # Ohms
L = 25e-3 # Henrios
C = 10e-6 # Faradios
I = 10 # Amperios -- Corriente inicial del inductor
V = 0 # Voltios -- Voltaje inicial del capacitor
tf = .004 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2


     Sistema Subamoritiguado:
     
R = 10 # Ohms
L = 1e-3 # Henrios
C = 35e-6 # Faradios
I = 0 # Amperios -- Corriente inicial del inductor
V = -75e-3 # Voltios -- Voltaje inicial del capacitor
tf = 1.5e-3 # Tiempo final(segundos)
S = 2 # Circuto en Paralelo = 1, Serie = 2
"""