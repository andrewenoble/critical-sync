# Python standard library modules.
import sys
from datetime import datetime

# Python package modules.
import scipy as sp
import cPickle as pk

# Local modules.
import libe as lb

# Global parameters and constants for local modules.
import global_params as gp
import global_consts as gc


#-------------------------- Define constants -------------------------------

R_MULT = gc.R_MULT


#-------------------------- Define parameters -------------------------------

# For directories.
load_dir = '../z_files/'
save_dir = 'z_files_pp/'

# MC sampling.
num_stats = int(1e5);
num_b = 100

len_b = num_stats / num_b
Gtij_list = [[0,0], [0,1], [1,1], [0,2], [1,2], [2,2], [0,3]]

# For writing files.
pickle_protocol = gp.pickle_protocol


#--------------------------------- Run  -------------------------------------

def main(dLR):
    
    d, L, R = dLR
    d = int(d)
    L = int(L)
    R /= R_MULT
    print d, L, R
    
    N = L * L 
    
    suffix = str(d) + "d_" + str(L) + "_" + str(int(round(R * R_MULT))) 
    
    # Load MC output.
    phi_array = abs(sp.loadtxt(load_dir + "phi_" + suffix + ".txt"))
    e_array = sp.loadtxt(load_dir + "e_" + suffix + ".txt")
    Gt00_array = sp.loadtxt(load_dir + "Gt00_" + suffix + ".txt")
    Gt01_array = sp.loadtxt(load_dir + "Gt01_" + suffix + ".txt")
    Gt11_array = sp.loadtxt(load_dir + "Gt11_" + suffix + ".txt")
    Gt02_array = sp.loadtxt(load_dir + "Gt02_" + suffix + ".txt")
    Gt12_array = sp.loadtxt(load_dir + "Gt12_" + suffix + ".txt")
    Gt22_array = sp.loadtxt(load_dir + "Gt22_" + suffix + ".txt")
    Gt03_array = sp.loadtxt(load_dir + "Gt03_" + suffix + ".txt") 
    
    Gt_arrays = [Gt00_array, Gt01_array, Gt11_array, Gt02_array, Gt12_array, 
        Gt22_array, Gt03_array]
          
    # Calcualate statistics.  
    
    ### Mean values.
    phi = sp.mean(phi_array)
    phi_var = sp.var(phi_array)
    phi_kurt = lb.binder(phi_array) 
    e = sp.mean(e_array)
    e_var = sp.var(e_array)  
    
    inverse_Gt_list = [1 / sp.mean(array) for array in Gt_arrays]
                            
    ### Autocorrelation times.
    tau_phi, tau_phi_err = lb.tau(phi_array)  
    tau_e, tau_e_err = lb.tau(e_array) 
    
    tau_Gt_list = [lb.tau(array)[0] for array in Gt_arrays]
    
    ### Error estimates based on autocorrelation times.
    phi_err = lb.err_ac(phi_array, tau_phi)
    e_err = lb.err_ac(e_array, tau_e)
    err_inverse_Gt_list = [lb.err_ac(array, tau) / sp.mean(array)**2
        for array, tau in zip(Gt_arrays, tau_Gt_list)]
          
    ### Error estimates based on blocking.
    phi_var_err = sp.std([sp.var(phi_array[i * len_b : (i+1) * len_b]) 
        for i in xrange(num_b)]) / sp.sqrt(num_b - 1)
    phi_kurt_err = sp.std([lb.binder(phi_array[i * len_b : (i+1) * len_b]) 
        for i in xrange(num_b)]) / sp.sqrt(num_b - 1)        
    e_var_err = sp.std([sp.var(e_array[i * len_b : (i+1) * len_b]) 
        for i in xrange(num_b)]) / sp.sqrt(num_b - 1) 
    
    ### Estimate xi and eta. 
    x = []  
    for i, j in Gtij_list:            
        x.append(lb.kfunc(L, i, j)) 
        
    x = sp.array(x)
        
    y = sp.array(inverse_Gt_list)   
    y0 = y[0]   
    y0_p =  1 / phi_var / N / N

    err = sp.array(err_inverse_Gt_list)
    err0 = err[0]   
    err0_p = phi_var_err / phi_var**2 / N / N 
    
    xi, xi_err, slope, intercept, dof, chi_sq_val = \
        lb.xi(N, x, y, err, y0, err0, y0_p, err0_p)
        
    # Record the results	            
    stats_list = [R, 
        phi, phi_err, 
        phi_var, phi_var_err, 
        phi_kurt, phi_kurt_err, 
        e, e_err, 
        e_var, e_var_err, 
        x, y, err, 
        xi, xi_err,
        slope, intercept,
        dof, chi_sq_val, 
        tau_phi, tau_phi_err, 
        sp.mean(Gt00_array) - (N * phi)**2]  

    # Save to file                    
    f = open(save_dir + 'stats_' + suffix, 'w') 
    pk.dump(stats_list, f, pickle_protocol)
    f.close()
    
    
if __name__ == "__main__":
    
    startTime = datetime.now()
        
    i = int(sys.argv[1])
    
    DLR_list = sp.loadtxt(load_dir + "DLR.txt")
    if i == 0: sp.savetxt(save_dir + "DLR.txt", DLR_list)

    main(DLR_list[i])
    
    print datetime.now() - startTime
