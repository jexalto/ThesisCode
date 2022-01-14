import scipy.special


def bezier(x, P):
    """Calculates values using a 1D Bezier curve
    
    Args:
        x (float/lst/numpy.ndarray): Points to be evaluated, must be between 0
            and 1
        P (lst/numpy.ndarray): Values of the control points
        
    Returns:
        s (lst/numpy.ndarray): Value on the Bezier curve
    """
    
    N = len(P)-1
    s = 0
    #Calculates values using explicit definition of a Bezier curve
    for i in range(N+1):
        s = s+scipy.special.binom(N, i)*(1-x)**(N-i)*x**i*P[i]
    
    return s
