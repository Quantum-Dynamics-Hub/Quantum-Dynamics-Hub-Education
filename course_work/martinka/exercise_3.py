def exp_taylor(arg):                                                       
    """
    exp(x) using Taylor series
    """
    result = 0
    for i in range(arg['n']+1):
        fact = 1
        for j in range(i):
            fact *= j+1
        result += pow(arg['x'],i)/fact
    return result

arg_dict = {'x':10, 'n':100}
exp_taylor(arg_dict)
