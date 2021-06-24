#!/user/bin/python
import math
def taylor(x, n):
    """
    this function computes exp(x) function using Taylor series
    
    -------------------------
    Arguments:
        x:int
        n:int
           the n value for taylor series
                      
    Returns:
        i: int
            the calculated value
          
    """
    i = 0
    a = 0
    while a <= n:
        i += x**a/math.factorial(a)
        a +=1
    return (i) 


print (taylor(3, 10))
