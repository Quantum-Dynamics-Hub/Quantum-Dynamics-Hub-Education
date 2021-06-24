import math

def myexp(x,n):
    e=1
    f=1
    for i in range(n):
        f*=i+1
        e+=x**(i+1)/f
    return e 

a=5
n=100
print('Python math library exp(',a,') = ',math.exp(a))
print('My power function exp(',a,') = ',myexp(a,n))
print('Difference is ', math.exp(a)-myexp(a,n))
