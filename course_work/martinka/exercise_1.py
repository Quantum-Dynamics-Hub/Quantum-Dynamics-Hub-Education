def pow_for(arg):
    """
    Ax^n, where n>0
    """
    if(arg['n']<0):
        raise ValueError('n = %i; n < 0' % (arg['n']))
    j = 1
    for i in range(arg['n']):
        j*=arg['x']
    result = j*arg['A']
    return result


arg_dict = {'A':1, 'x':3, 'n':-1}
pow_for(arg_dict)
