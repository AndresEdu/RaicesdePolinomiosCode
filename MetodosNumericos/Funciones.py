
def horner(grado, coeficientes, x):
    polinomio = coeficientes[grado]
    k = grado - 1
    while (k >= 0):
        polinomio = coeficientes[k] + (polinomio*x)
        k = k - 1
	# Al término de este ciclo WHILE, la variable polinomio tiene el valor del P(x)
    if polinomio==0:
        return polinomio+0.001
    else:
        return polinomio

def deriv(grado,coeficiente):
    k = 0
    derivados=[]
    while k <= grado:
        derivados.append((grado-k)*(coeficiente[k]))
        k += 1
    return derivados

def Newton(coeficientes,grado,tolerancia,xi,ernp,derivada):
    while (abs(ernp) >= tolerancia):
        numerador = horner(grado,coeficientes,xi)
        denominador = horner(grado,derivada,xi)
    # Obtener el valor actual x_(i+1)
        try: 
            x_imas1 = xi-(numerador/denominador)
        except:
            x_imas1 = xi-numerador

        # Cálculo del error: 
        ernp = ((x_imas1-xi)/x_imas1) * 100
        # El valor que se acaba de obtener es la i-ésima aproximación actual. Una vez obtenida, se toma como el valor xi y se vuelve a meter al polinomio. 
        xi = x_imas1
# Al término de este ciclo WHILE, la variable x_imas1 contiene la décima aproximación a la raíz del polinomio.

    return x_imas1

