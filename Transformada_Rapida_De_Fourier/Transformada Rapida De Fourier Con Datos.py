import numpy as np
import matplotlib.pyplot as plt

#Toma los datos de la memoria SD
# -------------------------------------------------------
archivo = "D:\SECCIONE.TXT"
f = open(archivo)
d = []
for one_line in f.readlines():
     d.append(list(map(float,one_line.strip().split())))
f.close()

y = [d[i][0] for i in range(len(d)-1)]  # Lista de las datos en la SD
t = [d[len(d)-1][0]/1000]               # Tiempo que se tardo en guardar los datos el Arduino
#-------------------------------------------------------

# Acomodo de datos para sacar la Transformas Rapida de Fourier (fft)
#-------------------------------------------------------
n = len(d)-1   # Numero de datos
t = t[0]       # Tiempo
dt = t/n       # Intervalos de tiempo dt = 1/fs = t/n
fs = 1/dt      # Frecuencia de muestreo fs = 1/dt

print(f"\nNumero de datos       (n) = {n}")
print(f"Tiempo                (t) = {t}")
print(f"Intervalo de tiempo  (∆t) = {dt}")
print(f"Frecuenca de muestro con ∆t (fs) = {fs}")

x = np.linspace(0,t,n) # Vector de tiempo
y = np.array(y)
#-------------------------------------------------------

# Ajuste de datos de un rango [0v,5v] a [-2.5v,2.5v]
#-------------------------------------------------------
s = 0
for i in range(len(y)):
     s += y[i]
p = s/n
y = y-p
#-------------------------------------------------------

# Calculo de la Transformada Rapida de Fourier (fft)
#-------------------------------------------------------
xt = np.arange(0,n,1)         # Vector para graficar fft sin frecuencia
xtn = xt                      
xtm = np.arange(0,n/2+1,1)    # Vector para graficar fft sin frecuencia mitad
xtmn = xtm

yt = np.fft.fft(y)       # Calculo de fft con numpy
ytr = np.abs(yt)         # Obtencian de la parte real de fft
yta = (2*ytr)/n          # Ajuste para tener la amplitud

ytm = np.fft.rfft(y)     # Calculo de fft con numpy mitad
ytrm = np.abs(ytm)       # Obtencian de la parte real de fft mitad
ytam = (2*ytrm)/n        # Ajuste para tener la amplitud mitad

#ytl = np.log((np.abs(yt)))    # Calcula el log_10 de la ttf
#ytl = 20*ytl

f = fs*xt/n                        # Calculo de la frecuecia
fm = fs*xtm/n                      # Calculo de la frecuecia mitad
fn = np.fft.fftfreq(len(x),dt)     # Calculo de la frecuencia con numpy
fnm = np.fft.rfftfreq(len(y),dt)   # Calculo de la frecuencia con numpy mitad
#-------------------------------------------------------

# Impresion de puntos maximos de la ttf
#-------------------------------------------------------
max_ytam = np.amax(ytam)
for i in range(len(ytam)):
     if ytam[i] == max_ytam:
          frecuencia = fm[i]
          max_xtm = xtm[i]
          
print(f"\nFrecuencia = {frecuencia}")
print(f"Max en x = {max_xtm}")

print(f"\nFrecuenca de muestro = {(n/2*max_xtm)/frecuencia}")
print(f"Frecuenca de N (fN) = {fs/2}")
#-------------------------------------------------------

# Grafica de la funcion muestreda
#-------------------------------------------------------
plt.style.use("fast") 
fig, g = plt.subplots(dpi = 110, figsize = (12.3,6.1))#, facecolor = "blue")
g.set_title(f"Señal\nf = {round(frecuencia,3)}", fontsize = 20)
g.set_xlabel("Tiempo (s)")
g.set_ylabel("x(t)")
g.grid(color = "blue")
g.axhline(color = "red")
g.axvline(color = "red")
g.plot(x,y,".-", color = "black")
g.plot(x,y,".", color = "blue")#,lw = 1, ms = 1)
#-------------------------------------------------------

# Grafica de la ttf completa
#-------------------------------------------------------
fig, g1 = plt.subplots(dpi = 110, figsize = (12.3,6.1))#, facecolor = "blue")
g1.set_title(f"Transformada Rapida de Fourier\nf = {round(frecuencia,3)}", fontsize = 20)
g1.set_ylabel("Amplitud")
g1.grid(color = "blue")
g1.axhline(color = "red")
g1.axvline(color = "red")
g1.plot(xtn,yta,"-", color = "black")
g1.plot(xtn,yta,".", color = "blue", label = "FFT")
#-------------------------------------------------------

# Grafica de la ttf completa
#-------------------------------------------------------
fig, g1 = plt.subplots(dpi = 110, figsize = (12.3,6.1))#, facecolor = "blue")
g1.set_title(f"Transformada Rapida de Fourier\nf = {round(frecuencia,3)}", fontsize = 20)
g1.set_ylabel("Amplitud")
g1.grid(color = "blue")
g1.axhline(color = "red")
g1.axvline(color = "red")
g1.plot(fn,yta,"-", color = "black")
g1.plot(fn,yta,".", color = "blue", label = "FFT")
#-------------------------------------------------------

# Grafica de la ttf la mitad de [0,pi]
#-------------------------------------------------------
#fig, g2 = plt.subplots(dpi = 110, figsize = (12.3,6.1))#, facecolor = "blue")
#g2.set_title(f"Transformada Rapida de Fourier\nf = {round(frecuencia,3)}", fontsize = 20)
#g2.set_ylabel("Amplitud")
#g2.grid(color = "blue")
#g2.axhline(color = "red")
#g2.axvline(color = "red")
#g2.plot((xtmn/n)*2*np.pi,ytam,"-", color = "black")
#g2.plot((xtmn/n)*2*np.pi,ytam,".", color = "blue", label = "FFT")
#-------------------------------------------------------

# Grafica de ttf con log_10 de [0,pi]
#-------------------------------------------------------
#fig, g3 = plt.subplots(dpi = 110, figsize = (12.3,6.1))#, facecolor = "blue")
#g3.set_title(f"Transformada Rapida de Fourier\nf = {round(frecuencia,3)}\n20*log_10()", fontsize = 15)
#g3.grid(color = "blue")
#g3.axhline(color = "red")
#g3.axvline(color = "red")
#g3.plot((xtn[1:n//2]/n)*2*np.pi,ytl[1:n//2],"-", color = "black")
#g3.plot((xtn[1:n//2]/n)*2*np.pi,ytl[1:n//2],".", color = "blue", label = "FFT")
#-------------------------------------------------------

plt.show()