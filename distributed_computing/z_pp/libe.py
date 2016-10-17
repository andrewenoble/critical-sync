# Python package modules.
import scipy as sp
import scipy.fftpack as fp
from scipy.optimize import curve_fit

# Global constants for local modules.
import global_consts as gc


#-------------------------- Define constants --------------------------------

K_TO_TAU_RATIO = gc.K_TO_TAU_RATIO

XI_NUM_FITTED_PARAMS = 2

#-------------------------- Define functions --------------------------------

xi_fit_func = lambda x, b, m: b + m * x

xi_fit_func_p = lambda x, a, b, c: a * (1 + (b * sp.sqrt(x))**c)

def kfunc(L, i, j):
    one = ((2 * sp.pi / L) * i)**2
    two = ((2 * sp.pi / L) * j)**2
    return one + two
    
    
def ac_fft(array):
    
    size = len(array)
    f = fp.fft(array - sp.mean(array), n=2 * size)
    ac = sp.real(fp.ifft(abs(f) ** 2))[:size]
    ac *= float(size) / (size - sp.arange(size))
    return ac / ac[0]
    

def tau(array):
    
    ac = ac_fft(array)
    tau = 1 / 2.
    for k, ack in enumerate(ac[1:]):
        tau += ack
        if ((k + 1) >= K_TO_TAU_RATIO * tau):
            break

    return tau, tau * sp.sqrt((2 * (2 * k + 1)) / float(len(array)))
    
    
def err_ac(array, tau):
    return sp.std(array) / sp.sqrt((len(array) - 1) / (2 * tau))


def xi_fit(f, x, y, err):

    popt, pcov = curve_fit(f, x, y, sigma=err)

    intercept = popt[0]
    slope = popt[1]
    
    e = intercept * sp.ones(len(x)) + slope * x
    chi_sq_val = sp.sum((y - e)**2 / err**2)
    
    dof = len(x) - XI_NUM_FITTED_PARAMS
    
    pcov = pcov * dof / chi_sq_val # See the explanation at http://www.physics.utoronto.ca/~phy326/python/curve_fit_to_data.py

    intercept_var = pcov[0, 0]
    slope_var = pcov[1, 1]
    cov = pcov[0, 1] 
    
    xi = sp.sqrt(slope / intercept)
    xi_err = (1 / sp.sqrt(4 * slope * intercept**3)) * \
        sp.sqrt(slope**2 * intercept_var + intercept**2 * slope_var
                - 2 * slope * intercept * cov)
                    
    return xi, xi_err, slope, intercept, dof, chi_sq_val

        
def xi_fit_pp(x, y, err):

    slope = (y[1] - y[0]) / (x[1] - x[0])
    intercept = (x[1] * y[0] - x[0] * y[1]) / (x[1] - x[0])
    
    xi = sp.sqrt((y[1] - y[0]) / (x[1] * y[0] - x[0] * y[1]))
    xi_err = (1 / (2. * xi)) * \
        sp.sqrt((x[0] - x[1])**2 * (err[1]**2 * y[0]**2 + err[0]**2 * y[1]**2) 
            / (x[1] * y[0] - x[0] * y[1])**4)
    '''
    slope_array = sp.array(
       [(y[1] - y[0] + err[1] + err[0]) / (x[1] - x[0]), 
        (y[1] - y[0] - err[1] - err[0]) / (x[1] - x[0])])
     
    intercept_array = sp.array(
        [y[0] - err[0] - x[0] * slope_array[0], 
         y[0] + err[0] - x[0] * slope_array[0]])
        
    xi_list = [sp.sqrt(slope / intercept) 
        for slope, intercept in zip(slope_array, intercept_array)] 
        
    xi = sp.mean(xi_list)       
    xi_err = sp.std(xi_list) #/ sp.sqrt(2 - 1)

    slope, intercept = sp.mean(slope_array), sp.mean(intercept_array)
    '''
    dof, chi_sq_val = 0, 0       
    
    return xi, xi_err, slope, intercept, dof, chi_sq_val
    
    
def xi(N, x, y, err, y0, err0, y0_p, err0_p):
    
    # Setup.
    upper = 6
    num_pts_list = [7, 6, 5, 4, 3]
    chi_sq_per_dof_95 = {
        1: 3.84 / 1, 2: 5.99 / 2, 3: 7.82 / 3, 4: 9.49 / 4, 5: 11.07 / 5}

    xi_index = 0
    xi_err_index = 1
    dof_index = 4
    chi_sq_val_index = 5
    
    # Assuming the disordered phase . . . 
    y[0] = y0
    err[0] = err0
    
    # Try the fit.                       
    results = xi_fit(xi_fit_func, x[:upper], y[:upper], err[:upper])
    chi_sq_per_dof = results[chi_sq_val_index] / results[dof_index]
    
    # Assuming the ordered phase . . .
    y[0] = y0_p
    err[0] = err0_p
    
    # Try the fit.
    results = xi_fit(xi_fit_func, x[:upper], y[:upper], err[:upper])
    chi_sq_per_dof_p = results[chi_sq_val_index] / results[dof_index]    
                    
                    
    # Let the two chi_sq_per_dof's calculated above determine whether 
    # we are in the ordered or disordered phase.                
    if chi_sq_per_dof < chi_sq_per_dof_p:
        
        # For the disordered phase:
        y[0] = y0
        err[0] = err0 
        
        lower = 0
        
        # Reduce the number of 1/Gt points till we have a good fit.                                  
        for upper in num_pts_list:
            
            # Try the fit.
            results = xi_fit(xi_fit_func, x[lower:upper], 
                y[lower:upper], err[lower:upper])
            chi_sq_per_dof = results[chi_sq_val_index] / results[dof_index]                        
    
            if chi_sq_per_dof < chi_sq_per_dof_95[results[dof_index]]:                 
                return results
                
        # We arrive here if none of the above fits are good.
        # In that case, simply estimate xi from the first two 1/Gt points.        
        results_2pt = xi_fit_pp(x[lower:lower+2], y[lower:lower+2], 
            err[lower:lower+2])  
        return results_2pt    
        
    else:
        
        # For the ordered phase:
        # let points the second and third 1/Gt points determine 
        # the value of xi. 
         
        lower = 1

        results_2pt = xi_fit_pp(x[lower:lower+2], y[lower:lower+2], 
            err[lower:lower+2])  
        return results_2pt   
                                            

# For low and high u behavior to match Ising model, m has to tend to
# 0 for high T, i.e. the asymptotic behavior of Binder's cumulant is
# different for high and low T only because the mean at high T is 0
# and the mean at low u is nonzero.
def binder(array):
    
    return sp.mean((array)**4)/sp.mean((array)**2)**2 - 3
