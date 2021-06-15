#!/user/bin/python
def my_pow(x, params):
    """
    this function computes A*x^n , for non-negative n values
    
    -------------------------
    Arguments:
        x:int
        params: dictionary
                the value of A and n in (A*x^n) should be given as a dictionary:
                for example: {'A':2, 'n':3}
                      
    Returns:
        i: int
            the calculated value
        
    if negative amount of n is given or if the input dictionary is not correct
    it will print an error message.    
    """
    # first check if the doctionary is correct
    
    if params.keys() == {'A', 'n'}:
        A = params['A']
        n = params ['n']
        if n<0:
            return ('Error message: negative value of n is not accepted')
        else:
            i = 1
            for k in range (n):
                i = i * x
            i = A * i
            return (i)
    else:
        return ('Error message: the input disctionary is not correct')

        
params = {'A':4}
print (my_pow(5, params))

params = {'A':4, 'm':1}
print (my_pow(5, params))

params = {'A':4, 'n':-1}
print (my_pow(5, params))

params = {'A':4, 'n':0}
print (my_pow(5, params))

params = {'A':4, 'n':2}
print (my_pow(5, params))


