import math

def finint(func, arg):
    """
    Solution for trivial cases
    function to compute a definite integral of an aribtrary function using Simpson's 1/3 rule
    """
    arg_a, arg_b, arg_ab = arg.copy(), arg.copy(), arg.copy()
    del arg_a['ab']
    del arg_b['ab']
    del arg_ab['ab']
    arg_a['x'] = arg['ab'][0]
    arg_b['x'] = arg['ab'][1]
    arg_ab['x'] = (arg['ab'][0]+arg['ab'][1])/2
    res = (arg_b['x']-arg_a['x'])*(func(arg_a)+4*func(arg_ab)+func(arg_b))/6
    return res

def arb_func1(var): #function a*x^2
    res = var['a']*var['x']**2
    return res

def arb_func2(var): #function a*exp(b*x)
    res = var['a']*math.exp(var['b']*var['x'])
    return res

def arb_func3(var): #function d*x^3+c*x^2-b*x+a
    res = var['d']*var['x']**3+var['c']*var['x']**2-var['b']*var['x']+var['a']
    return res

arg1 = {'ab':[10,20], 'x':2, 'a':10}
arg2 = {'ab':[1,2], 'x':2, 'a':10, 'b':2}
arg3 = {'ab':[10,20], 'x':2, 'a':10, 'b':10, 'c':10, 'd':10}

print(finint(arb_func1,arg1))
print(finint(arb_func2,arg2))
print(finint(arb_func3,arg3))
