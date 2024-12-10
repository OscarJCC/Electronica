import numpy as np

def trunc(cf,num):

     if type(num) == float or type(num) == int:
          n = str(round(float(num),cf))
          l = 0
          nu =""
          for j in n:
               if j == "e":
                    nu += j
                    l = 1   
               elif l == 1:
                    nu += j
                    
          if (len(n)-len(nu)) <= cf:
               num = float(str(n[:len(n)-len(nu)]))
          elif num < 0 and num - int(num) == 0:
               num = float(str(n[:cf+1]))
          elif num < 0 and num - int(num) != 0:
               num = float(str(n[:cf+2]))
          elif num - int(num) != 0:
               num = float(str(n[:cf+1]))
          elif num - int(num) == 0:
               num = float(str(n[:cf]))

          for j in str(num):
               if j == "e":
                    l = 0
                    break
                    
          if l == 1: 
               num = str(num)
               num += nu
               num = float(num)
               l = 0
          
     elif type(num) == complex:
          rnum = round(float(num.real),cf)
          rnumi = round(float(num.imag),cf)
          rni = str(rnum)
          ri = str(rnumi)
                    
          if rnum < 0 and rnum - int(rnum) == 0:
               rnum = float(str(rni[:cf+1]))
          elif rnum < 0 and rnum - int(rnum) != 0:
               rnum = float(str(rni[:cf+2]))
          elif rnum - int(rnum) != 0:
               rnum = float(str(rni[:cf+1]))
          elif rnum - int(rnum) == 0:
               rnum = float(str(rni[:cf]))
          
          if rnumi < 0 and rnumi - int(rnumi) == 0:
               rnumi = float(str(ri[:cf+1]))
          elif rnumi < 0 and rnumi - int(rnumi) != 0:
               rnumi = float(str(ri[:cf+2]))
          elif rnumi - int(rnumi) != 0:
               rnumi = float(str(ri[:cf+1]))
          elif rnumi - int(rnumi) == 0:
               rnumi = float(str(ri[:cf]))
          
          num = rnum + (rnumi)*1j
     
     elif type(num) == list:
          nume = []
          for i in range(len(num)):
               n = str(round(float(num[i]),cf))
               
               l = 0
               nu =""
               for j in n:
                    if j == "e":
                         nu += j
                         l = 1   
                    elif l == 1:
                         nu += j
               
               if num[i] < 0 and num[i] - int(num[i]) == 0:
                    nume.append(float(str(n[:cf+1])))
               elif num[i] < 0 and num[i] - int(num[i]) != 0:
                    nume.append(float(str(n[:cf+2])))
               elif num[i] - int(num[i]) != 0:
                    nume.append(float(str(n[:cf+1])))
               elif num[i] - int(num[i]) == 0:
                    nume.append(float(str(n[:cf])))
                    
               for j in str(nume[i]):
                    if j == "e":
                         l = 0
                    
               if l == 1: 
                    nume[i] = str(nume[i])
                    nume[i] += nu
                    nume[i] = float(nume[i])
                    l = 0
          
     if type(num) == list:
          return nume   
     else:       
          return num

def PrintMat(name,matriz,c,f,cf):
     if cf > 5:
          print(f"\n{name}:\n")
          for i in range(c):
               for j in range(f):
                    if j == f-1:
                         print("{:^3}{:^24}{:^3}".format(" ¦ ",trunc(cf,matriz[i][j])," ¦ "),end=" ")
                    else:
                         print("{:^3}{:^24}".format(" ¦ ",trunc(cf,matriz[i][j])),end=" ")
               print("\n")
     elif cf > 10:
          print(f"\n{name}:\n")
          for i in range(c):
               for j in range(f):
                    if j == f-1:
                         print("{:^3}{:^24}{:^3}".format(" ¦ ",trunc(cf,matriz[i][j])," ¦ "),end=" ")
                    else:
                         print("{:^3}{:^24}".format(" ¦ ",trunc(cf,matriz[i][j])),end=" ")
               print("\n")
     else:
          print(f"\n{name}:\n")
          for i in range(c):
               for j in range(f):
                    if j == f-1:
                         print("{:^3}{:^24}{:^3}".format(" ¦ ",trunc(cf,matriz[i][j])," ¦ "),end=" ")
                    else:
                         print("{:^3}{:^24}".format(" ¦ ",trunc(cf,matriz[i][j])),end=" ")
               print("\n")


v = 15

r1 = 10e3
r2 = 10e2
r3 = 10e3
r4 = 47e2
r5 = 56e1
r6 = 10e4

I = 0

lr = [r1,r2,r3,r4,r5,r6]

ec_sn = [-((1/r1)+(1/r5)),1/r1,1/r5,-1/r3]

ec_a = [1/r1,-((1/r1)+(1/r2)+(1/r6)),1/r2,0]

ec_c = [1/r5,1/r2,-((1/r2)+(1/r4)+(1/r5)),0]

ec_aux = [1,0,0,-1]

mat = [ec_sn,ec_a,ec_c,ec_aux]

sol = [0,0,0,15]

A = np.array(mat)
B = np.array(sol)

Ainv = np.linalg.inv(A)

solf = np.dot(Ainv,B)

lv = list(solf)
     
va = solf[0]
vb = solf[1]
vc = solf[2]
vd = solf[3]

ir1 = (va-vb)/r1
ir2 = (vb-vc)/r2
ir3 = -vd/r3
ir4 = -vc/r4
ir5 = (va-vc)/r5
ir6 = -vb/r6

Iv = -(ir1+ir5) 

li = [ir1,ir2,ir3,ir4,ir5,ir6]

PrintMat("Matriz",mat,4,4,24)

print("Solucion:\n","\t"," ¦ ",sol[0]," ¦ ",sol[1]," ¦ ",sol[2]," ¦ ",sol[3]," ¦ ","\n")

for i in range(len(lr)):
     print(f"R{i+1} = ",lr[i])

print("\n")

print("Va = ",va)
print("Vb = ",vb)
print("Vc = ",vc)
print("Vd = ",vd)

print("\n")

for i in range(len(li)):
     print(f"IR{i+1} = ",li[i])
     I += li[i]

print("\n")

print("It = ",I)

print("\n")

print("Iv =",Iv)

print("Suma de Corriente It + Iv = ",I+Iv)

