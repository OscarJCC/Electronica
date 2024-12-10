import numpy as np
import matplotlib.pyplot as plt

n = 2**10 # Numero de datos
fs = 2**9 # Frecuencia de muestreo fs = 1/dt
dt = 1/fs # Intervalos de tiempo dt = 1/fs = t/n
t = n*dt  # Tiempo

f1 = 5    # Frecuencia 1
f2 = 12.5 # Frecuencia 2
f3 = 20   # Frecuencia 3
f4 = 155  # Frecuencia 4

A1 = 2    # Amplitud 1
A2 = 1    # Amplitud 2
A3 = 1.5  # Amplitud 3
A4 = 0.5  # Amplitud 4

x = np.linspace(0,t,n) # Vector de tiempo
y = A1*np.sin(f1*2*np.pi*x)+A2*np.sin(f2*2*np.pi*x)+A3*np.sin(f3*2*np.pi*x)+A4*np.sin(f4*2*np.pi*x)
#y = A3*np.sin(f3*2*np.pi*x)

xt = np.arange(0,n,1)         # Vector para graficar fft sin frecuencia
xtm = np.arange(0,n/2+1,1)    # Vector para graficar fft sin frecuencia mitad

yt = np.fft.fft(y)  # Calculo de fft con numpy
ytr = np.abs(yt)    # Obtencian de la parte real de fft
yta = (2*ytr)/n     # Ajuste para tener la amplitud

ytm = np.fft.rfft(y)  # Calculo de fft con numpy mitad
ytrm = np.abs(ytm)    # Obtencian de la parte real de fft mitad
ytam = (2*ytrm)/n     # Ajuste para tener la amplitud mitad

ytl = np.log(np.abs(yt))
ytl = 20*ytl

f = fs*xt/n                        # Calculo de la frecuecia
fm = fs*xtm/n                      # Calculo de la frecuecia mitad
fn = np.fft.fftfreq(len(x),dt)     # Calculo de la frecuencia con numpy
fnm = np.fft.rfftfreq(len(y),dt)   # Calculo de la frecuencia con numpy mitad

# Grafica de funcion
plt.style.use("dark_background") 
fig, g = plt.subplots(dpi = 110, figsize = (12.3,6.1), facecolor = "blue")
g.set_title(f"Funcion", fontsize = 20)
g.set_xlabel("Tiempo (s)")
g.set_ylabel("x(t)")
g.grid(color = "blue")
g.axhline(color = "red")
g.axvline(color = "red")
g.plot(x,y, color = "yellow")
#g.legend(edgecolor = "red", fontsize = 15)

# Grafica de Amplitud vs vector xt
fig, g1 = plt.subplots(dpi = 110, figsize = (12.3,6.1), facecolor = "blue")
g1.set_title("Transformada Rapida de Fourier\nEspectro de Amplitud", fontsize = 20)
g1.set_ylabel("Amplitud")
g1.grid(color = "blue")
g1.axhline(color = "red")
g1.axvline(color = "red")
g1.plot(xt,yta, color = "yellow", label = "FFT")
#g1.legend(edgecolor = "red", fontsize = 15)

# Grafica de Amplitud vs vector xt mitad
fig, g2 = plt.subplots(dpi = 110, figsize = (12.3,6.1), facecolor = "blue")
g2.set_title("Transformada Rapida de Fourier\nEspectro de Amplitud", fontsize = 20)
g2.set_ylabel("Amplitud")
g2.grid(color = "blue")
g2.axhline(color = "red")
g2.axvline(color = "red")
g2.plot(xtm,ytam, color = "yellow", label = "FFT")
#g2.legend(edgecolor = "red", fontsize = 15)

# Grafica de Amplitud vs Frecuencia
fig, g3 = plt.subplots(dpi = 110, figsize = (12.2,6.1), facecolor = "blue")
g3.set_title("Transformada Rapida de Fourier\nAmplitud vs Frecuencia", fontsize = 20)
g3.set_xlabel("Frecuencia")
g3.set_ylabel("Amplitud")
g3.grid(color = "blue")
g3.axhline(color = "red")
g3.axvline(color = "red")
g3.plot(f,yta, color = "yellow", label = "FFT")
#g3.legend(edgecolor = "red", fontsize = 15)

# Grafica de Amplitud vs Frecuencia
fig, g4 = plt.subplots(dpi = 110, figsize = (12.2,6.1), facecolor = "blue")
g4.set_title("Transformada Rapida de Fourier\nAmplitud vs Frecuencia", fontsize = 20)
g4.set_xlabel("Frecuencia")
g4.set_ylabel("Amplitud")
g4.grid(color = "blue")
g4.axhline(color = "red")
g4.axvline(color = "red")
g4.plot(fm,ytam, color = "yellow", label = "FFT")
#g4.legend(edgecolor = "red", fontsize = 15)

fig, g5 = plt.subplots(dpi = 110, figsize = (12.3,6.1), facecolor = "blue")
g5.set_title("Transformada Rapida de Fourier\n20*log_10()", fontsize = 20)
g5.grid(color = "blue")
g5.axhline(color = "red")
g5.axvline(color = "red")
g5.plot(xt[1:n],ytl[1:n], color = "yellow", label = "FFT")
#g5.legend(edgecolor = "red", fontsize = 15)

fig, g5 = plt.subplots(dpi = 110, figsize = (12.3,6.1), facecolor = "blue")
g5.set_title("Transformada Rapida de Fourier\n20*log_10() con frecuencia", fontsize = 20)
g5.grid(color = "blue")
g5.axhline(color = "red")
g5.axvline(color = "red")
g5.plot(f[1:n],ytl[1:n], color = "yellow", label = "FFT")
#g5.legend(edgecolor = "red", fontsize = 15)

plt.show()
