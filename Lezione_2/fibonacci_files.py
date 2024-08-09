import numpy as np
def fibonacci(n):
#FIBONACCI  Fibonacci sequence
#   f = FIBONACCI(n) generates the first n Fibonacci numbers.

     f = np.zeros((n,))
     f[0] = 1.0
     f[1] = 1.0
     for k in range(2,n):
         f[k] = f[k-1] + f[k-2]
     fn=f[-1]   
     return fn


def fib_ric(n):
   
    if n <= 2:   
        return 1
    else:
       f = fib_ric(n-1) + fib_ric(n-2)
       return f

