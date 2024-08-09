from sign import sign
import math
import numpy as np

def bisezione(fname, a, b, tolx):
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
 if sign(fa)*sign(fb)>=0:
     print("Non è possibile applicare il metodo di bisezione \n")
     return None, None,None

 it = 0
 v_xk = []

 
 maxit = math.ceil(math.log2((b - a) / tolx))-1

 
 while abs(b - a) > tolx:
    xk = a+(b-a)/2
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)
    if fxk==0:
      return xk, it, v_xk

    if sign(fa)*sign(fxk)>0:  #continua su [xk,b]
      a = xk
      fa=fxk
    elif sign(fxk)*sign(fb)>0:   #continua su [a,xk]
      b = xk
      fb=fxk

 
 return xk, it, v_xk


def corde(fname,m,x0,tolx,tolf,nmax):
    
     # m è il coefficiente angolare della retta che rimane fisso per tutte le iterazioni
        xk=[]
        fx0=fname(x0)
        d=fx0/m
        x1=x0-d
        fx1=fname(x1)
        xk.append(x1)
        it=1
        
        while it<nmax and  abs(fx1)>=tolf and abs(d)>=tolx*abs(x1) :
           x0=x1
           fx0=fname(x0)
           d=fx0/m
           '''
           #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
           '''
           x1=x0-d  
           fx1=fname(x1)
           it=it+1
         
           xk.append(x1)
          
        if it==nmax:
            print('raggiunto massimo numero di iterazioni \n')
            
        
        return x1,it,xk



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
 fa=fname(a);
 fb=fname(b);
 if sign(fa)*sign(fb)>=0:
     print("Non è possibile applicare il metodo di falsa posizione \n")
     return None, None,None

 it = 0
 v_xk = []
 
 fxk=10

 
 while it < maxit and abs(b - a) > tolx and abs(fxk) > tolf:
    xk = a-fa*(b-a)/(fb-fa)
    v_xk.append(xk)
    it += 1
    fxk=fname(xk)
    if fxk==0:
      return xk, it, v_xk

    if sign(fa)*sign(fxk)>0:  #continua su [xk,b]
      a = xk
      fa=fxk
    elif sign(fxk)*sign(fb)>0:   #continua su [a,xk]
      b = xk
      fb=fxk

 
 return xk, it, v_xk



def newton_mod(fname,fpname,m,x0,tolx,tolf,nmax):
  #m: molteplicità della radice
        xk=[]
        fx0=fname(x0)
        if abs(fpname(x0))<=np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
            print(" derivata prima nulla in x0")
            return None, None,None

        d=fx0/fpname(x0)
        x1=x0-m*d
        
        fx1=fname(x1)
        xk.append(x1)
        it=1
        
        while it<nmax and  abs(fx1)>=tolf and abs(d)>=tolx*abs(x1) :
           x0=x1
           fx0=fname(x0)
           if abs(fpname(x0))<=np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
                print(" derivata prima nulla in x0")
                return None, None,None
           d=fx0/fpname(x0)
           '''
           #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) ed è tangente alla funzione f(x) nel punto (xi.f(xi))  e l'asse x
           '''
           x1=x0-m*d  
           fx1=fname(x1)
           it=it+1
         
           xk.append(x1)
          
        if it==nmax:
            print('raggiunto massimo numero di iterazioni \n')
            
        
        return x1,it,xk



def newton(fname,fpname,x0,tolx,tolf,nmax):
  
        xk=[]
        fx0=fname(x0)
        if abs(fpname(x0))<=np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
            print(" derivata prima nulla in x0")
            return None, None,None
        
        d=fx0/fpname(x0)
        x1=x0-d
        
        fx1=fname(x1)
        xk.append(x1)
        it=1
        
        while it<nmax and  abs(fx1)>=tolf and abs(d)>=tolx*abs(x1) :
           x0=x1
           fx0=fname(x0)
           if abs(fpname(x0))<=np.spacing(1): #Se la derivata prima e' pià piccola della precisione di macchina stop
                print(" derivata prima nulla in x0")
                return None, None,None
           d=fx0/fpname(x0)
           '''
           #x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) ed è tangente alla funzione f(x) nel punto (xi.f(xi))  e l'asse x
           '''
           x1=x0-d  
           fx1=fname(x1)
           it=it+1
         
           xk.append(x1)
          
        if it==nmax:
            print('raggiunto massimo numero di iterazioni \n')
            
        
        return x1,it,xk

def secanti(fname,xm1,x0,tolx,tolf,nmax):
        xk=[]
        fxm1=fname(xm1)
        fx0=fname(x0)
        d=fx0*(x0-xm1)/(fx0-fxm1)
        x1=x0-d
        xk.append(x1)
        fx1=fname(x1)
        it=1
       
        while it<nmax and abs(fx1)>=tolf and abs(d)>=tolx*abs(x1):
            xm1=x0
            x0=x1
            fxm1=fname(xm1)
            fx0=fname(x0) 
            d=fx0*(x0-xm1)/(fx0-fxm1)
            x1=x0-d
            fx1=fname(x1)
            xk.append(x1)
            it=it+1
           
       
        if it==nmax:
           print('Secanti: raggiunto massimo numero di iterazioni \n')
        
        return x1,it,xk