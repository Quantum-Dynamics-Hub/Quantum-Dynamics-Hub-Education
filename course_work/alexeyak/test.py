print("Hi")

def my_pow1(x, params):
    """
    Computes x^n , for non-negative n values 
    """
    A, B, D, b = 10.0, 0.1, 3.0, 3
    
    A = params["A"]
    B = params["B"]
    D = params["D"]
    n = params["n"]

    print( type(n) )
    res = A*x**n + B + D*math.log(x)

    return res 

