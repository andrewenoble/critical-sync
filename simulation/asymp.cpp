#include "nr3.h" 
#include "ran.h"
#include "deviates.h"
    
    
/*--------------------------- Define constants ------------------------------*/

const Int PREC = 12;
const Doub TOL = 1.0e-12; // A local variable will be set to 0 if less than TOL.
const Int CHAR_SIZE = 1000; 


/*---------------------------- Setup output ---------------------------------*/

// Open an output stream.
ofstream f;
        
// Specify a directory for saving files.
char save_dir[CHAR_SIZE] = "../simulation_output/\0"; // C strings must be null terminated.
    
// Create a string that will contain save_dir + filename.
char buf[CHAR_SIZE] = {0};


/*--------------------- Define simulation parameters ------------------------*/

// These parameters govern the Monte Carlo sampling of the asymptotic 
// synchronization order parameter configurations.
Int num_sweeps_burnin = int(0); // No burnin needed if "x_init.txt" is loaded.
Int num_sweeps = int(2e4); // Should be an EVEN number.
Int num_samples = int(2e2); // Should be an EVEN number.

// These parameters define the stochastic coupled map lattice model that we
// want to simulate.
Doub lam = 0.15;
Doub eps = 0.20;
Doub R = 2.166135;
Int L = 256; // L defines the linear dimension of a 2D square network.


/*------------------- Define derived simulation parameters ------------------*/

Int N = L * L; // N is the total number of nodes in the 2D square network.
Int pause_between_samples = num_sweeps / num_samples;
Doub coupling = eps / 4.; // Coupling is symmetric among the four nearest
                          // neighbors on the 2D square network.


/*------------------------------ Utilities ------------------------------*/
    
// Define a function that will save an array of floats to file.  
void save_to_file(VecDoub &X, Char *name) {
    
    f.open(name);
    
    for (Int i = 0; i < X.size(); ++i) f << setprecision(PREC) << X[i] << "\n"; 
    
    f.close();
    
}


// Define a function to setup the nearest-neighbor look-up table, "nn".
void neighbor(MatInt &nn, Int &L, Int &N){  
    	
	for (int k = 0; k < N ; k++) {
	    
		nn[k][0] = (k / L) * L + (k + 1) % L;  //+x
		nn[k][1] = (k / L) * L + (k - 1 + L) % L;  //-x		
		nn[k][2] = (k + L) % N; //+y
		nn[k][3] = (k - L + N) % N; //-y

	}
	
}


// Linear diffusive iteractions of each node on the network with its 
// nearest-neighbors is simulated using a nearest-neighbor look-up table, "nn".
void ax(MatInt &nn, VecDoub &x, VecDoub &y, Doub &coupling, Int &N) {

    for (Int i = 0; i < N; ++i) {
                                
        y[i] = (1. - 4. * coupling) * x[i] + coupling * (x[ nn[i][0] ] 
            + x[ nn[i][1] ] + x[ nn[i][2] ] + x[ nn[i][3] ]);
        
    }
  
}


// Define two functions that iterate the stochastic map for each node in the
// network.  The first function overwrites the vector "x".
void map(VecDoub &x, Doub &lam, Doub &R, Int &N, SNormaldev &SNormaldev) {
    
    Doub value;
        
    for (Int i = 0; i < N; ++i) {
        
        value = x[i];
        
        if (value > TOL) {
            
            value *= exp(R * (1.0 - value));
            value = (value > TOL) ? value : 0.0;
            value += SNormaldev.dev() * lam * value;
            x[i] = (value > TOL) ? value : 0.0;
        
        } else {
            
            x[i] = 0.0;
        
        }
    }
}

// The second function does not overwrite "x".
void map(VecDoub &x, VecDoub &y, Doub &lam, Doub &R, Int &N, 
    SNormaldev &SNormaldev) {
  
    Doub value;
    
    for (Int i = 0; i < N; ++i) {
        
        value = x[i];
        
        if (value > TOL) {
            
            value *= exp(R * (1.0 - value));
            value = (value > TOL) ? value : 0.0;
            value += SNormaldev.dev() * lam * value;
            y[i] = (value > TOL) ? value : 0.0;
        
        } else {
            
            y[i] = 0.0;
        
        }
    }
}


/*------------------------------- Run --------------------------------*/

// Run a Monte Carlo simluation of a stochastic coupled map lattice.
Int mc() {
        
    // Read in an initial network configuration obtained from early simulations.
    // (Almost any initial condition could be chosen in principle, but the user
    // would need to set "num_sweeps_burnin" to a large number and run for a
    // a long time in order to observe the critical behavior that emerges 
    // in the asymptotic dynamics.)
    ifstream in("x_init.txt");
    
    VecDoub x(N, 0.0);
    
    for (Int i = 0; i < N; i++) {
        
        in >> x[i];
        
    }
            
    // Generate the nearest neighbor table.
    MatInt nn(N, 4, 0);
    neighbor(nn, L, N);
    
    // I fix the seed here to a particular value.  
    // Any choice of a positive definite integer seed will generate the 
    // desired critical dynamics.
    SNormaldev SNormaldev(12345);

    // Iterate the stochastic process to burn in.
    Int i; 
    VecDoub y(N, 0.0);
    
    for (i = 0; i < num_sweeps_burnin; ++i) {

        // Iterate twice per "sweep".
        map(x, lam, R, N, SNormaldev);
        ax(nn, x, y, coupling, N);
        
        map(y, lam, R, N, SNormaldev);
        ax(nn, y, x, coupling, N);
                
    }
    
    // Iterate in the asymptotic regime to gather snapshots of the 
    // synchronization order parameter configuration.
    VecDoub dummy(N, 0.0), m(N, 0.0);
    Int j, k = -1;       
                        
    for (i = 0; i < num_sweeps; ++i) {
        
        // Iterate twice per "sweep".
        map(x, lam, R, N, SNormaldev);
        ax(nn, x, y, coupling, N);
        
        map(y, dummy, lam, R, N, SNormaldev);
        ax(nn, dummy, x, coupling, N);
        
        // Pause regularly to calculate and save the synchronization order 
        // parameter configuration.
        if ((i % pause_between_samples) == 0) { 
            
            for (j = 0; j < N; ++j) m[j] = (x[j] - y[j]) / 2;
                    
            sprintf(buf, "%sm_%d.txt", save_dir, ++k);
            save_to_file(m, buf);  
            
        }                            
                    
    }
                      
    return 1;
    
}


int main(int argc, char** argv) {
                
    // Run that Monte Carlo simulation.                
    Int return_value = mc();
    
    return 0;
    
}