def pow_for(arg):
    """
    Ax^n, where n>0
    """
    if(arg['n']<0):
        raise ValueError('n = %i; n < 0' % (arg['n']))
    j,i = 1,1
    while(i<=arg['n']):
        j*=arg['x']
        i+=1
    result = j*arg['A']
    return result

arg_dict = {'A':1, 'x':3, 'n':4}
pow_for(arg_dict)
