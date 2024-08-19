import math
import numpy as np
import numpy.linalg as npl
import scipy.linalg as spl
import sympy as sym
from SolveTriangular import *


def metodo_bisezione(fname, a, b, tolx,tolf):
 """
 Implementa il metodo di bisezione per il calcolo degli zeri di un'equazione non lineare.

 Parametri:
  f: La funzione da cui si vuole calcolare lo zero.
  a: L'estremo sinistro dell'intervallo di ricerca.
  b: L'estremo destro dell'intervallo di ricerca.
  tol: La tolleranza di errore.

 Restituisce:
  Lo zero approssimato della funzione, il numero di iterazioni e la lista di valori intermedi.
 """

 fa=fname(a)
 fb=fname(b)
 if math.copysign(1,fa)*math.copysign(1,fb)>0 : 
     print("Non è possibile applicare il metodo di bisezione \n")
     return None, None,None

 it = 0
 v_xk = []

 maxit = math.ceil(math.log((b - a) / tolx) / math.log(2))-1
 
 
 while it <= maxit and abs(b-a) >= tolx: 
    xk = a + (b-a)/2
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)


    # Controlla se fxk è sufficientemente vicino a zero
    if abs(fxk) < tolf:
      return xk, it, v_xk

    # Se xk è a sinistra dello 0, allora muovo a a destra 
    if math.copysign(1,fa)*math.copysign(1,fxk)>0:   
      a = xk
      fa=fxk
    # Se xk è a destra dello 0, allora muovo b a sinistra
    elif math.copysign(1,fxk)*math.copysign(1,fb)>0:    
      b = xk
      fb=fxk

 return xk, it, v_xk


def falsi(fname, a, b, maxit, tolx,tolf):
 """
 Implementa il metodo di falsa posizione per il calcolo degli zeri di un'equazione non lineare.

 Parametri:
  f: La funzione da cui si vuole calcolare lo zero.
  a: L'estremo sinistro dell'intervallo di ricerca.
  b: L'estremo destro dell'intervallo di ricerca.
  tol: La tolleranza di errore.

 Restituisce:
  Lo zero approssimato della funzione, il numero di iterazioni e la lista di valori intermedi.
 """
 fa=fname(a)
 fb=fname(b)
 
 if  math.copysign(1,fa)*math.copysign(1,fb)>=0 : 
      print("Non è possibile applicare il metodo di falsa posizione \n")
      return None, None,None

 it = 0
 v_xk = []
 
 fxk=10

 
 while it <= maxit and abs(b-a)>tolx and abs(fxk) > tolf: 
    xk = a - fa*(b-a)/(fb -fa)
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)
    if abs(fxk) <= tolf:
      return xk, it, v_xk

    # 
    if math.copysign(1,fa)*math.copysign(1,fxk)>0:
        a = xk
        fa = fxk
    elif math.copysign(1,fxk)*math.copysign(1,fb)>0:    
        b = xk
        fb = fxk
    

 
 return xk, it, v_xk


def corde(fname,m,x0,tolx,tolf,nmax):
  """
  Implementa il metodo delle corde per il calcolo degli zeri di un'equazione non lineare.

  Parametri:
  fname: La funzione da cui si vuole calcolare lo zero.
  m: coefficiente angolare della retta che rimane fisso per tutte le iterazioni
  tolx: La tolleranza di errore tra due iterati successivi
  tolf: tolleranza sul valore della funzione
  nmax: numero massimo di iterazione

  Restituisce:
  Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
  """
  xk=[]
  fx0=fname(x0)
  d= fx0/m
  x1=x0 -d
  fx1=fname(x1)
  xk.append(x1)
  it=1
      
  while it<nmax and  abs(fx1)>=tolf and abs(d)>=tolx*abs(x1) :
    x0 = x1
    fx0= fname(x0)
    d= fx0/m
    '''
    #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
    (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
    '''
    x1 = x0 -d
    fx1=fname(x1)
    xk.append(x1)
    it=it+1

  return x1,it,xk

def newton(fname,fpname,x0,tolx,tolf,nmax):
  """
  Implementa il metodo di Newton per il calcolo degli zeri di un'equazione non lineare.

  Parametri:
    fname: La funzione di cui si vuole calcolare lo zero.
    fpname: La derivata prima della funzione di  cui si vuole calcolare lo zero.
    x0: iterato iniziale
    tolx: La tolleranza di errore tra due iterati successivi
    tolf: tolleranza sul valore della funzione
    nmax: numero massimo di iterazione

  Restituisce:
    Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
  """ 
  xk=[]
  fx0=fname(x0)
  if abs(fpname(x0)) <= np.spacing(1): 
      print(" derivata prima nulla in x0")
      return None, None,None

  d=fx0/fpname(x0)
  x1= x0-d

  fx1=fname(x1)
  xk.append(x1)
  it=1

  while abs(fx1-fx0) >= tolf and it  < nmax and abs(d) >= tolx * abs(x1):
      x0= x1
      fx0= fname(x0)
      if abs(fpname(x0))<=np.spacing(1) :
          print(" derivata prima nulla in x0")
          return None, None,None
      

      d=fx0/fpname(x0)
      
      x1=x0 - d
      fx1=fname(x1)
      it=it+1
    
      xk.append(x1)
    
  if it==nmax:
      print('raggiunto massimo numero di iterazioni \n')
      

  return x1,it,xk

def secanti(fname,xm1,x0,tolx,tolf,nmax):
  """
  Implementa il metodo delle secanti per il calcolo degli zeri di un'equazione non lineare.

  Parametri:
    fname: La funzione di cui si vuole calcolare lo zero.
    xm1, x0: primi due iterati
    tolx: La tolleranza di errore tra due iterati successivi
    tolf: tolleranza sul valore della funzione
    nmax: numero massimo di iterazione

  Restituisce:
    Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
  """
  xk=[]
  fxm1=fname(xm1)
  fx0=fname(x0)
  d= fx0 * (x0-xm1)/(fx0-fxm1)
  x1= x0 - d
  xk.append(x1)
  fx1=fname(x1)
  it=1
  
  while it<nmax and abs(fx1)>=tolf and abs(d)>=tolx*abs(x1):
      xm1= x0
      x0 = x1
      fxm1= fname(xm1)
      fx0=fname(x0)
      d=fx0 * (x0-xm1)/(fx0-fxm1)
      x1=x0 - d
      fx1=fname(x1)
      xk.append(x1)
      it=it+1

  return x1,it,xk

def newton_mod(fname,fpname,m,x0,tolx,tolf,nmax):
  """
  Implementa il metodo di Newton modificato da utilizzato per il calcolo degli zeri di un'equazione non lineare
  nel caso di zeri multipli.

 Parametri:
  fname: La funzione di cui si vuole calcolare lo zero.
  fpname: La derivata prima della funzione di  cui si vuole calcolare lo zero.
   m: molteplicità della radice
  x0: iterato iniziale
  tolx: La tolleranza di errore tra due iterati successivi
  tolf: tolleranza sul valore della funzione
  nmax: numero massimo di iterazione

 Restituisce:
  Lo zero approssimato della funzione, il numero di iterazioni e la lista degli iterati intermedi.
  """ 
 
  
  xk=[]
  fx0=fname(x0)
  if abs(fpname(x0)) <= np.spacing(1): 
      print(" derivata prima nulla in x0")
      return None, None,None

  d=fx0/fpname(x0)
  x1= x0-m*d

  fx1=fname(x1)
  xk.append(x1)
  it=1

  while abs(fx1) >= tolf and it  < nmax and abs(d) >= tolx * abs(x1):
      x0= x1
      fx0= fname(x0)
      if abs(fpname(x0))<=np.spacing(1) :
          print(" derivata prima nulla in x0")
          return None, None,None
      

      d=fx0/fpname(x0)
      
      x1=x0 - m*d
      fx1=fname(x1)
      it=it+1
    
      xk.append(x1)
    
  if it==nmax:
      print('raggiunto massimo numero di iterazioni \n')
      

  return x1,it,xk
    
def stima_ordine(xk,iterazioni):
  #Vedi dispensa allegata per la spiegazione

  k=iterazioni-4
  p=np.log(abs(xk[k+2]-xk[k+3])/abs(xk[k+1]-xk[k+2]))/np.log(abs(xk[k+1]-xk[k+2])/abs(xk[k]-xk[k+1]));
  
  ordine=p
  return ordine

def my_newtonSys(fun, jac, x0, tolx, tolf, nmax):

  """
  Funzione per la risoluzione del sistema F(x)=0
  mediante il metodo di Newton.

  Parametri
  ----------
  fun : funzione vettoriale contenente ciascuna equazione non lineare del sistema.
  jac : funzione che calcola la matrice Jacobiana della funzione vettoriale.
  x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
  tolx : float
    Parametro di tolleranza per l'errore assoluto.
  tolf : float
    Parametro di tolleranza per l'errore relativo.
  nmax : int
    Numero massimo di iterazioni.

  Restituisce
  -------
  x : array
    Vettore soluzione del sistema (o equazione) non lineare.
  it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
  Xm : array
    Vettore contenente la norma dell'errore relativo tra due iterati successivi.
  """

  matjac = jac(x0)
  if np.linalg.det(matjac) == 0:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None, None,None

  s = - npl.solve(matjac, fun(x0))
  # Aggiornamento della soluzione
  it = 1
  x1 = x0 + s
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

  while it <= nmax and npl.norm(fx1, 1) >= tolf and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1):
    x0 = x1
    it += 1
    matjac = jac(x0)
    if npl.det(matjac) == 0:
            print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
            return None, None,None

   
    s = - npl.solve(matjac, fun(x0))

    # Aggiornamento della soluzione
    x1 = x0 + s
    fx1 = fun(x1)
    Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

  return x1, it, Xm


def my_newtonSys_corde(fun, jac, x0, tolx, tolf, nmax):

  """
  Funzione per la risoluzione del sistema f(x)=0
  mediante il metodo di Newton.

  Parametri
  ----------
  fun : stringa
    Nome del file contenente la funzione non lineare.
  jac : stringa
    Nome del file contenente la matrice Jacobiana della funzione.
  x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
  tolx : float
    Parametro di tolleranza per l'errore assoluto.
  tolf : float
    Parametro di tolleranza per l'errore relativo.
  nmax : int
    Numero massimo di iterazioni.

  Restituisce
  -------
  x : array
    Vettore soluzione del sistema (o equazione) non lineare.
  it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
  Xm : array
      Vettore contenente la norma dell'errore relativo tra due iterati successivi.
  """

  matjac = jac(x0)  #Utilizzo per tutte le iterazioni la matrice Jacobiana valutata nell'ierato iniziale, senza mai aggiornarla
  if np.linalg.det(matjac) == 0:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None, None,None
  s = -np.linalg.solve(matjac, fun(x0))
  # Aggiornamento della soluzione
  it = 1
  x1 = x0 + s
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]

  while it <= nmax and np.linalg.norm(fx1, 1) >= tolf and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1):
    x0 = x1
    it += 1
   
   
    if np.linalg.det(matjac) == 0:
        print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
        return None, None,None
    
    # Risolvo il sistema lineare avente come matrice dei coefficienti la
    # matrice Jacobiana e come termine noto la Funzione vettoriale F valutata
    # in x0
    
    s = -np.linalg.solve(matjac, fun(x0))

    # Aggiornamento della soluzione
    x1 = x0 + s
    fx1 = fun(x1)
    Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

  return x1, it, Xm

def my_newtonSys_sham(fun, jac, x0, tolx, tolf, nmax):

  """
  Funzione per la risoluzione del sistema f(x)=0
  mediante il metodo di Newton.

  Parametri
  ----------
  fun : stringa
    Nome del file contenente la funzione non lineare.
  jac : stringa
    Nome del file contenente la matrice Jacobiana della funzione.
  x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
  tolx : float
    Parametro di tolleranza per l'errore assoluto.
  tolf : float
    Parametro di tolleranza per l'errore relativo.
  nmax : int
    Numero massimo di iterazioni.

  Restituisce
  -------
  x : array
    Vettore soluzione del sistema (o equazione) non lineare.
  it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
  Xm : array
      Vettore contenente la norma dell'errore relativo tra due iterati successivi.
  """

  matjac = jac(x0)
  if np.linalg.det(matjac) == 0:
    print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
    return None,None,None

  s = -np.linalg.solve(matjac, fun(x0))
  # Aggiornamento della soluzione
  it = 1
  x1 = x0 + s
  fx1 = fun(x1)

  Xm = [np.linalg.norm(s, 1)/np.linalg.norm(x1,1)]
  update=10  #Numero di iterazioni durante le quali non si aggiorna la valutazione dello Jacobiano nell'iterato attuale
  while it <= nmax and np.linalg.norm(fx1, 1) >= tolf and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1):
    x0 = x1
    it += 1
    if it%update==0:   #Valuto la matrice di iterazione nel nuovo iterato ogni "update" iterazioni
        matjac=jac(x0)
   
        if np.linalg.det(matjac) == 0:
           print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
           return None,None,None
        else:
         # Risolvo il sistema lineare avente come matrice dei coefficienti la
        # matrice Jacobiana valutatata nell'iterato aggiornato x0  e come termine noto la Funzione vettoriale F valutata
        # in x0
           s = -np.linalg.solve(matjac, fun(x0))
    else:
         # Risolvo il sistema lineare avente come matrice dei coefficienti la
        # matrice Jacobiana non aggiornata e come termine noto la Funzione vettoriale F valutata
        # in x0
           s = -np.linalg.solve(matjac, fun(x0))

    # Aggiornamento della soluzione
    x1 = x0 + s
    fx1 = fun(x1)
    Xm.append(np.linalg.norm(s, 1)/np.linalg.norm(x1,1))

  return x1, it, Xm


def my_newton_minimo(gradiente, Hess, x0, tolx, tolf, nmax):
    """
    DA UTILIZZARE NEL CASO IN CUI CALCOLATE DERIVATE PARZIALI PER GRADIENTE ED HESSIANO SENZA UTILIZZO DI SYMPY
    
    Funzione di Newton-Raphson per calcolare il minimo di una funzione in più variabili.
    
    Parametri
    ----------
    gradiente : 
        Nome della funzione che calcola il gradiente della funzione non lineare.
    Hess :  
        Nome della funzione che calcola la matrice Hessiana della funzione non lineare.
    x0 : array
        Vettore contenente l'approssimazione iniziale della soluzione.
    tolx : float
        Parametro di tolleranza per l'errore assoluto.
    tolf : float
        Parametro di tolleranza per l'errore relativo.
    nmax : int
        Numero massimo di iterazioni.
    
    Restituisce
    -------
    x : array
        Vettore soluzione del sistema (o equazione) non lineare.
    it : int
        Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
    Xm : array
        Vettore contenente la norma del passo ad ogni iterazione.
    """
    
    # Calcolo dell'Hessiana nel vettore iniziale x0
    matHess = Hess(x0)
    
    # Controllo se la matrice Hessiana è invertibile
    if np.linalg.det(matHess) == 0:
        print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
        return None, None, None
    
    # Calcolo il gradiente nel vettore iniziale x0
    grad_fx0 = gradiente(x0)
    
    # Calcolo del passo iniziale s risolvendo il sistema lineare Hess * s = -gradiente
    s = -np.linalg.solve(matHess, gradiente(x0))
    
    # Aggiornamento della soluzione iniziale con il nuovo passo
    it = 1  # Inizializzo il contatore delle iterazioni
    x1 = x0 + s  # Aggiorno la soluzione
    grad_fx1 = gradiente(x1)  # Calcolo il gradiente nel nuovo punto x1
    Xm = [np.linalg.norm(s, 1)]  # Salvo la norma del passo iniziale in una lista
    
    # Ciclo iterativo per aggiornare la soluzione fino alla convergenza o al massimo numero di iterazioni
    while it <= nmax and np.linalg.norm(grad_fx1, 1) >= tolf and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1):
        
        # Aggiornamento della soluzione corrente
        x0 = x1
        it += 1
        
        # Calcolo dell'Hessiana e del gradiente nel nuovo punto x0
        matHess = Hess(x0)
        grad_fx0 = grad_fx1
        
        # Controllo se la matrice Hessiana è invertibile
        if np.linalg.det(matHess) == 0:
            print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
            return None, None, None
        
        # Calcolo del nuovo passo s risolvendo Hess * s = -gradiente
        s = -np.linalg.solve(matHess, grad_fx0)
        
        # Aggiornamento della soluzione con il nuovo passo
        x1 = x0 + s
        
        # Calcolo del gradiente nel nuovo punto x1
        grad_fx1 = gradiente(x1)
        
        # Stampa della norma del passo corrente (utile per il debug o per monitorare la convergenza)
        print(np.linalg.norm(s, 1))
        
        # Aggiungo la norma del passo corrente alla lista Xm
        Xm.append(np.linalg.norm(s, 1))
    
    # Ritorno la soluzione finale, il numero di iterazioni e la lista delle norme dei passi
    return x1, it, Xm

def newton_minimo_MOD(gradiente, Hess, x0, tolx, tolf, nmax):

  """
  Funzione di newton-raphson per calcolare il minimo di una funzione in più variabili, modificato nel caso in cui si utilizzando sympy 
  per calcolare Gradinete ed Hessiano. Rispetto alla precedente versione cambia esclusivamente il modo di valutare il vettore gradiente e la matrice Hessiana in un punto 
  Parametri
   ----------
  fun : 
    Nome della funzione che calcola il gradiente della funzione non lineare.
  Hess :  
    Nome della funzione che calcola la matrice Hessiana della funzione non lineare.
  x0 : array
    Vettore contenente l'approssimazione iniziale della soluzione.
  tolx : float
    Parametro di tolleranza per l'errore assoluto.
  tolf : float
    Parametro di tolleranza per l'errore relativo.
  nmax : int
    Numero massimo di iterazioni.

  Restituisce
  -------
  x : array
    Vettore soluzione del sistema (o equazione) non lineare.
  it : int
    Numero di iterazioni fatte per ottenere l'approssimazione desiderata.
  Xm : array
    Vettore contenente la norma del passo ad ogni iterazione.
  """

    
  matHess = np.array([[Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                      [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])]])
 

  gradiente_x0=np.array([gradiente[0](x0[0], x0[1]),gradiente[1](x0[0], x0[1])])
   
  if np.linalg.det(matHess) == 0:
    print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
    return None, None, None
      
  s = -np.linalg.solve(matHess, gradiente_x0)
  # Aggiornamento della soluzione
  it = 1
  x1 = x0 + s
  grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
  Xm = [np.linalg.norm(s, 1)]
  
  while it <= nmax and np.linalg.norm(grad_fx1, 1) >= tolf and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1):
     
    x0 = x1
    it += 1
    matHess = np.array([[Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                      [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])]])
    grad_fx0=grad_fx1
      
    if np.linalg.det(matHess) == 0:
       
      print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
      return None, None, None
      

    # Risolvo il sistema lineare avente come matrice dei coefficienti la
    # matrice Hessiana e come termine il vettore gradiente calcolato nell'iterato precedente
    # in x0
                #NB: in fx1 è memorizzato il gradiente nell'iterato attuale
    s = -np.linalg.solve(matHess, grad_fx0)
     
    # Aggiornamento della soluzione
    x1 = x0 + s
    #Aggiorno il gradiente per la prossima iterazione 
    grad_fx1=np.array([gradiente[0](x1[0],x1[1]),gradiente[1](x1[0],x1[1])])
    Xm.append(np.linalg.norm(s, 1))

  return x1, it, Xm

# A = M - N
# A = D + E + F
# con D diagonale
# E triangolare inferiore con diagonale = 0, 
# F triangolare superiore con diagonale = 0
#
# quindi M = D
#        N = -(E+F)
#
def prof_jacobi(A,b,x0,toll,it_max):
    errore=1000
    D=np.diag(A)
    n=A.shape[0]
    M_inv=np.diag(1/D)
    E=np.tril(A,-1)
    F=np.triu(A,1)
    N=-(E+F)
    T=np.dot(M_inv,N)
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.max(np.abs(autovalori))
    print("raggio spettrale jacobi", raggiospettrale)
    it=0
    
    er_vet=[]
    while it < it_max and errore >= toll:
        x=(b+np.dot(N,x0))/D.reshape(n,1)
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

# A = M - N
# A = D + E + F
# con D diagonale
# E triangolare inferiore con diagonale = 0, 
# F triangolare superiore con diagonale = 0
#
# quindi M = D
#        N = -(E+F)
#
def jacobi(A,b,x0,toll,it_max):
    errore=1000
    d=np.diag(A)
    n=A.shape[0]
    invM=np.diag(1/d)
    E= np.tril(A, -1)
    F= np.triu(A, 1)
    N= - (E + F)
    T= invM@N
    autovalori=np.linalg.eigvals(T)
    raggiospettrale= np.max(np.abs(autovalori))
    print("raggio spettrale jacobi", raggiospettrale)
    it=0
    
    er_vet=[]
    while errore >= toll and it < it_max:
        x= T@x0 + invM@b
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet

# A = M - N
# A = D + E + F
# con D diagonale
# E triangolare inferiore con diagonale = 0, 
# F triangolare superiore con diagonale = 0
#
# quindi M = E+D
#        N = -F
#
def gauss_seidel(A,b,x0,toll,it_max):
    errore=1000
    D= np.diag(A)
    E=np.tril(A,-1)
    F=np.triu(A,1)
    M= D+E
    N= -F
    T= np.dot(npl.inv(M), N)
    autovalori=np.linalg.eigvals(T)
    raggiospettrale= np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel ",raggiospettrale)
    it=0
    er_vet=[]
    while it < it_max and errore >= toll:
        temp= b - F@x0
        x, flag = Lsolve(M, temp)
        errore=np.linalg.norm(x-x0)/np.linalg.norm(x)
        er_vet.append(errore)
        x0=x.copy()
        it=it+1
    return x,it,er_vet


def LUsolve(P,L,U,b):
    pb=np.dot(P,b)
    y,flag=Lsolve(L,pb)
    if flag == 0:
         x,flag=Usolve(U,y)
    else:
        return [],flag

    return x,flag

def gauss_seidel_sor(A,b,x0,toll,it_max,omega):
    errore=1000
    d=np.diag(A)
    D=np.diag(d)
    Dinv= np.diag(1/d)
    E= np.tril(A,-1)
    F= np.triu(A,1)
    M=D+E
    N=-F
    Momega=D+omega*E
    Nomega=(1-omega)*D-omega*F
    T=np.dot(np.linalg.inv(Momega),Nomega)
    autovalori=np.linalg.eigvals(T)
    raggiospettrale=np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)
    
    
    it=0
    xold=x0.copy()
    xnew=x0.copy()
    er_vet=[]
    while it < it_max and errore >= toll : 
        temp= b-F@xold
        xtilde,flag=Lsolve(M,temp)
        xnew=(1-omega)*xold+omega*xtilde
        errore=np.linalg.norm(xnew-xold)/np.linalg.norm(xnew)
        er_vet.append(errore)
        xold=xnew.copy()
        it=it+1
    return xnew,it,er_vet

def steepestdescent(A,b,x0,itmax,tol):
 
    n,m=A.shape
    if n!=m:
        print("Matrice non quadrata")
        return [],[]
    
    
   # inizializzare le variabili necessarie
    x = x0

     
    r = A@x-b  # residuo
    p = -r     # opposto del residuo
    it = 0
  
    errore=np.linalg.norm(r)/np.linalg.norm(b)
    vec_sol=[]
    vec_sol.append(x)
    vet_r=[]
    vet_r.append(errore)
     
# utilizzare il metodo del gradiente per trovare la soluzione
    while errore>= tol and it< itmax:
        it=it+1
        Ap=A@p
       
        alpha = (r.T@r)/(p.T@Ap)
                
        x = x + alpha*p  #aggiornamento della soluzione nella direzione opposta a quella del gradiente: alpha mi dice dove fermarmi 
        #nella direzione del gradiente affinche F(xk+t p ) <F(xk)
        
         
        vec_sol.append(x)
        r=r+alpha*Ap
        errore=np.linalg.norm(r)/np.linalg.norm(b)
        vet_r.append(errore)
        p = -r #Direzione opposta alla direzione del gradiente
        
     
    return x,vet_r,vec_sol,it


def conjugate_gradient(A,b,x0,itmax,tol):
    n,m=A.shape
    if n!=m:
        print("Matrice non quadrata")
        return [],[]
    
    
   # inizializzare le variabili necessarie
    x = x0
    
    r = A@x-b
    p = -r
    it = 0
    nb=np.linalg.norm(b)
    errore=np.linalg.norm(r)/nb
    vec_sol=[]
    vec_sol.append(x0)
    vet_r=[]
    vet_r.append(errore)
# utilizzare il metodo del gradiente coniugato per calcolare la soluzione
    while errore >= tol and it< itmax:
        it=it+1
        Ap=A@p
        alpha = -(r.T@p)/(Ap.T@p)
        x = x + alpha *p
        vec_sol.append(x)
        rtr_old=r.T@r
        r=r+alpha*Ap
        gamma=r.T@r/rtr_old
        errore=np.linalg.norm(r)/nb
        vet_r.append(errore)
        p = -r+gamma*p  #La nuova direzione appartiene al piano individuato da -r e p. gamma è scelto in maniera tale che la nuova direzione
        #sia coniugata rispetto alla direzione precedente( che geometricamente significa che punti verso il centro)
   
    
    return x,vet_r,vec_sol,it
    
def eqnorm(A,b):
 
    G=A.T@A
     
    f=A.T@b
    L=spl.cholesky(G,lower=True)
   
    y,flag= Lsolve(L,f)
    if flag==0:
        x,flag=Usolve(L.T,y)
    
    
    return x
def qrLS(A,b):
#Risolve un sistema sovradeterminato con il metodo QR-LS
    n=A.shape[1]  # numero di colonne di A
    Q,R=spl.qr(A)
    h= Q.T@b
    x,flag=Usolve(R[0:n,:],h[0:n])
    residuo=np.linalg.norm(h[n:])**2
    return x,residuo

def SVDLS(A,b):
    #Risolve un sistema sovradeterminato con il metodo SVD-LS
    m,n=A.shape  #numero di righe e  numero di colonne di A
    U,s,VT=spl.svd(A)  #Attenzione : Restituisce U, il numpy-array 1d che contiene la diagonale della matrice Sigma e VT=VTrasposta)
    #Quindi 
    V=VT.T
    thresh=np.spacing(1)*m*s[0] ##Calcolo del rango della matrice, numero dei valori singolari maggiori di una soglia
    k=np.count_nonzero(s>thresh)
    print("rango=",k)
    d= U.T@b
    d1= d[0:k].reshape(k,1)
    s1= s[:k].reshape(k,1)
    #Risolve il sistema diagonale di dimensione kxk avene come matrice dei coefficienti la matrice Sigma
    c= d1/s1
    x=V[:,:k]@c 
    residuo=np.linalg.norm(d[k:])**2
    return x,residuo

def plagr(xnodi,j):
  """
  Restituisce i coefficienti del j-esimo pol di
  Lagrange associato ai punti del vettore xnodi
  """
  xzeri=np.zeros_like(xnodi)
  n=xnodi.size
  if j==0:
      xzeri=xnodi[1:n]
  else:
      xzeri=np.append(xnodi[0:j], xnodi[j+1:n])
  
  num= np.poly(xzeri)
  den= np.polyval(num, xnodi[j])
  
  p=num/den
  
  return p

def InterpL(x, y, xx):
  """
  %funzione che determina in un insieme di punti il valore del polinomio
  %interpolante ottenuto dalla formula di Lagrange.
  % DATI INPUT
  %  x  vettore con i nodi dell'interpolazione
  %  f  vettore con i valori dei nodi 
  %  xx vettore con i punti in cui si vuole calcolare il polinomio
  % DATI OUTPUT
  %  y vettore contenente i valori assunti dal polinomio interpolante
  %
  """
  n=x.size
  m=xx.size
  L=np.zeros((m,n))
  
  for j in range(n):
    p= plagr(x,j)
    L[:,j]= np.polyval(p,xx)

  return L@y