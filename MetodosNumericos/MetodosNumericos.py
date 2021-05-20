
#Arciniega Arellano Andres Eduardo
#Programa que lea un polinomio de grado n > 0 y localice sus puntos notables.
import math
import sys
import cmath as cmath
from math import inf


#Algoritmo para sacar Bairstown
def Bairstow(coeficientes, r, s, grado, raiz, tolerancia):
	if(grado<1):
		return None
	if((grado==1) and (coeficientes[1]!=0)):
		raiz.append(float(-coeficientes[0])/float(coeficientes[1]))
		return None
	if(grado==2):
		D = (coeficientes[1]**2.0)-(4.0)*(coeficientes[2])*(coeficientes[0])
		if(D<0):
			X1 = (-coeficientes[1] - cmath.sqrt(D))/(2.0*coeficientes[2])
			X2 = (-coeficientes[1] + cmath.sqrt(D))/(2.0*coeficientes[2])
		else:
			X1 = (-coeficientes[1] - math.sqrt(D))/(2.0*coeficientes[2])
			X2 = (-coeficientes[1] + math.sqrt(D))/(2.0*coeficientes[2])
		raiz.append(X1)
		raiz.append(X2)
		return None
	n = len(coeficientes)
	b = [0]*len(coeficientes)
	c = [0]*len(coeficientes)
	b[n-1] = coeficientes[n-1]
	b[n-2] = coeficientes[n-2] + r*b[n-1]
	i = n - 3
	while(i>=0):
		b[i] = coeficientes[i] + r*b[i+1] + s*b[i+2]
		i = i - 1
	c[n-1] = b[n-1]
	c[n-2] = b[n-2] + r*c[n-1]
	i = n - 3
	while(i>=0):
		c[i] = b[i] + r*c[i+1] + s*c[i+2]
		i = i - 1
	Din = ((c[2]*c[2])-(c[3]*c[1]))**(-1.0)
	r = r + (Din)*((c[2])*(-b[1])+(-c[3])*(-b[0]))
	s = s + (Din)*((-c[1])*(-b[1])+(c[2])*(-b[0]))
	if(abs(b[0])> tolerancia or abs(b[1])>tolerancia):
		return Bairstow(coeficientes,r,s,grado,raiz, tolerancia)
	if (grado>=3):
		Dis = ((r)**(2.0))+((4.0)*(1.0)*(s))
		X1 = (r - (cmath.sqrt(Dis)))/(2.0)
		X2 = (r + (cmath.sqrt(Dis)))/(2.0)
		raiz.append(X1)
		raiz.append(X2)
		return Bairstow(b[2:],r,s,grado-2,raiz, tolerancia)

#Algoritmo de Horner
def horner(grado, coeficientes, x):
    polinomio = coeficientes[grado]
    k = grado - 1
    while (k >= 0):
        polinomio = coeficientes[k] + (polinomio*x)
        k = k - 1
    return polinomio

#Obtencion de la derivada
def deriv(Grado, COEF):
    k = 0
    derivado = []
    if((len(COEF)) == 1):
        derivado.append(0)
    else:
        while k < Grado:
            derivado.append((Grado - k)*(COEF[k]))
            k += 1
    return derivado

#--------------Inicio del Programa----------------

#Se ingresa el grado del Polinomio
Grado = int(input("Grado del polinomo: "))
Coeficientes = []

if Grado < 1 :
    print("Ingrese un grado valido ")
    sys.exit()

#Ingresamos los coeficientes de nuestro polinomio
aux = Grado;
for i in range(Grado+1):
    coeficiente = float(input("Ingresa el coeficiente x^"+str(aux)+": "))
    aux = aux-1
    Coeficientes.append(coeficiente)

#Ingresamos las cifras significativas
CifrasSignif = int(input("Cifras significativas: "))
print("")
print("")
#Obtenemos la primera derivada del polinomio
PrimeraDerivada = deriv(Grado, Coeficientes)         
#Obtenemos la segunda derivada del polinomio
SegundaDerivada = deriv(Grado-1, PrimeraDerivada)     

#Valores necesarios en el Algoritmo de Bairstow
r = Coeficientes[-1]/Coeficientes[0]
s = r
r1 = PrimeraDerivada[-1]/PrimeraDerivada[0]
r2 = SegundaDerivada[-1]/SegundaDerivada[0]
s1 = r1
s2 = r2
tolerancia = 0.5*pow(10, 2-CifrasSignif)   

Coeficientes.reverse()
PrimeraDerivada.reverse()
SegundaDerivada.reverse()

#Listas necesarias
RaicesPolinomio = []
RaicesPrimDeriv = []
RaicesSegDeriv = []

#Sacamos las raices imaginarios con el Algoritmo de Bairstow
Bairstow(Coeficientes, r, s, Grado, RaicesPolinomio, tolerancia)
Bairstow(PrimeraDerivada, r1, s1, Grado- 1, RaicesPrimDeriv, tolerancia)
Bairstow(SegundaDerivada, r2, s2, Grado- 2, RaicesSegDeriv, tolerancia)

#Separar raices complejas de las reales
RaicesPrimDeriv = [x for x in RaicesPrimDeriv if type(x) is not complex]
RaicesSegDeriv = [x for x in RaicesSegDeriv if type(x) is not complex]

# Máximos y mínimos
maximos = []    
minimos = []

c = 1/pow(10, CifrasSignif)

if(len(RaicesPrimDeriv) != 0):
    for x in RaicesPrimDeriv:
        primer = horner(Grado-1, PrimeraDerivada, (x-c))
        segundo = horner(Grado-1, PrimeraDerivada, (x+c))
    
    if ((primer < 0) and (segundo > 0)):
            minimos.append(complex(x, horner(Grado, Coeficientes, x)))
            minimos.append((x, horner(Grado, Coeficientes, x)))

    elif ((primer > 0) and (segundo < 0)):
            maximos.append(complex(x, horner(Grado, Coeficientes, x)))
            maximos.append((x, horner(Grado, Coeficientes, x)))
else:
    maximos.append("No hay Maximos")
    minimos.append("No hay Minimos")

#Puntosa de Inflexion
inflexion = []

if(len(RaicesSegDeriv) != 0) and Grado > 2:
        for x in RaicesSegDeriv:
            inflexion.append((x, horner(Grado, Coeficientes, x)))
else:
    inflexion.append("No hay puntos de inflexion")

#Puntos donde es creciente y decreciente
crece = []
decrece = []
RaicesPrimDeriv.sort()

if Grado == 1:
    crece.append("No crece")
    decrece.append("No decrece")
else:
    if(len(RaicesPrimDeriv) != 0):
        for i in range(1, len(RaicesPrimDeriv)):
            if horner(Grado-1, PrimeraDerivada, RaicesPrimDeriv[i]-c) < 0:
                decrece.append((RaicesPrimDeriv[i-1], RaicesPrimDeriv[i]))
            else:
        #crece.append(complex(RaicesPrimDeriv[i-1], RaicesPrimDeriv[i]))
                crece.append((RaicesPrimDeriv[i-1], RaicesPrimDeriv[i]))

        if horner(Grado-1, PrimeraDerivada, RaicesPrimDeriv[-1]+c) < 0:
            #decrece.append(complex(RaicesPrimDeriv[-1], inf))
             decrece.append((RaicesPrimDeriv[-1], inf))
        else:
            #crece.append(complex(RaicesPrimDeriv[-1], inf))
             crece.append((RaicesPrimDeriv[-1], inf))

        if horner(Grado-1, PrimeraDerivada, RaicesPrimDeriv[0] - c < 0):
             decrece.append((-inf, RaicesPrimDeriv[0]))
        else:
             crece.append((-inf, RaicesPrimDeriv[0]))
    else:
        crece.append("No crece")
        decrece.append("No decrece")

#intervalos de concavidad
arriba = []
abajo = []
RaicesSegDeriv.sort()

if(len(RaicesSegDeriv) != 0) and Grado > 2:
    for i in range(1, len(RaicesSegDeriv)):
        if horner(Grado-2, SegundaDerivada, RaicesSegDeriv[i]-c) < 0:
            abajo.append((RaicesSegDeriv[i-1], RaicesSegDeriv[i]))
        else:
            arriba.append((RaicesSegDeriv[i-1], RaicesSegDeriv[i]))

    if horner(Grado-2, SegundaDerivada, RaicesSegDeriv[-1] + c) < 0:
        abajo.append((RaicesSegDeriv[-1], inf))
    else:
        arriba.append((RaicesSegDeriv[-1], inf))

    if horner(Grado-2, SegundaDerivada, RaicesSegDeriv[0] - c) < 0:
        abajo.append((-inf, RaicesSegDeriv[0]))
    else:
        arriba.append((-inf, RaicesSegDeriv[0]))
else:
    arriba.append("No es concava hacia arriba \n")
    abajo.append("No es concava hacia abajo \n")

#Imprimimos en Pantalla
print("RAICES del polinomio:")
print("")

for i in RaicesPolinomio:
    print(i)

print("")
print("Minimos:")
print("")
for x in minimos:
    print(x, end=", ")

print("")
print("Maximos:")
print("")
for x in maximos:
    print(x, end=", ")

print("")
print("Puntos de Inflexion:")
print("")
for x in inflexion:
    print(x, end=", ")

print("")
print("Decreciente en: ")
print("")
for x in decrece:
    print(x, end=", ")

print("")
print("Creciente en: ")
print("")
for x in crece:
    print(x, end=", ")

print("")
print("Concava hacia abajo en: ")
print("")
for x in abajo:
    print(x, end=", ")

print("")
print("Concava hacia arriba en: ")
print("")
for x in arriba:
    print(x, end=", ")

