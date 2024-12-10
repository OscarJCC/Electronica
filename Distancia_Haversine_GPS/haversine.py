
import numpy as np
def cV(lx,ly): # Funcion -> Concatenacion Vertical
     l = []
     
     for i in range(len(lx)):
          l.append([lx[i],ly[i]])
          
     return np.array(l)

def haversine(nlat,nlong,vlat,vlong):
     R = 6371000
     phi1 = nlat*(np.pi/180.0);
     phi2 = vlat*(np.pi/180.0);
     delta1 = (vlat-nlat)*(np.pi/180.0);
     delta2 = (vlong-nlong)*(np.pi/180.0);
     a = pow(np.sin(delta1/2.0),2.0) + np.cos(phi1)*np.cos(phi2)*pow(np.sin(delta2/2.0),2.0);
     c = 2*np.arctan2(np.sqrt(a),np.sqrt(1.0-a));
     d = R*c;
     
     return d

lat = np.array([25.413740,25.413965,25.414155,25.413984,25.413724,25.413896])
long = np.array([-100.986915,-100.987396,-100.987335,-100.986915,-100.986915,-100.986648])
distancia = np.array([0.0]*len(lat))

for i in range(len(lat)):
     if i != 0:
          distancia[i] = haversine(lat[i],long[i],lat[i-1],long[i-1])

posicion = cV(lat,long)
print(posicion)
print(distancia)
