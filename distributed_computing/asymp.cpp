#include "nr3.h" 
#include "sparse.h"  
#include "ran.h"
#include "deviates.h"
#include "sort.h"
#include "fourier.h"
#include "fourier_ndim.h"
#include <time.h>


/*-------------------------- Start the clock ---------------------------*/

clock_t start = clock();
    
    
/*----------------------------- Constants ------------------------------*/

const Int PREC = 12;
const Doub TOL = 1.0e-12; // Floor for each local population.
const Int LAM_MULT = int(1e6); // Precision with which R-values are recorded.
const Int CHAR_SIZE = 1000; // String size.

/*--------------------------------- IO ---------------------------------*/

// Open an output stream.
ofstream f;
        
// Specify a directory for saved files.
char save_dir[CHAR_SIZE] = "z_files/\0"; // C strings must be null terminated.
    
// Create a string that will contain save_dir + filename.
char buf[CHAR_SIZE] = {0};


/*-------------------------- Model parameters ---------------------------*/

// MC sampling
Int num_sweeps_burnin = int(1e7);
Int num_sweeps = int(1e8);
Int num_x = int(1e1);
Int num_running_updates = int(1e6);
Int num_stats = int(1e5);

// Dynamics
Doub mu0 = 0.5;
Doub sig20 = 0.01;

Doub lam_low = 0.135;
Doub lam_high = 0.150;
Int lam_num = 72;

Doub eps = 0.10;
Doub R = 2.3;

Int L_low = 16;
Doub L_pow_inc = 2.0;
Int L_num = 3;

Int D_low = 2;    
Int D_inc = 1;    
Int D_num = 1;    


/*-------------------------- Derived parameters -------------------------*/

Int x_record = num_sweeps - num_x / 2;
Int running_stats_record = num_sweeps / num_running_updates;
Int list_stats_record = num_sweeps / num_stats;

/*------------------------------ Utilities ------------------------------*/

// Save arrays.
void save_to_file(MatDoub &X, Char *name) {
    
    f.open(name);
    
    for (Int i = 0; i < X.nrows(); ++i) {
        
        for (Int j = 0; j < X.ncols(); ++j) 
        
            f << setprecision(PREC) << X[i][j] << " "; 

	f << "\n";
	
    }
    
    f.close();
    
}
    

void save_to_file(VecDoub &X, Char *name) {
    
    f.open(name);
    
    for (Int i = 0; i < X.size(); ++i) f << setprecision(PREC) << X[i] << "\n"; 
    
    f.close();
    
}


void save_to_file(Uint &seed, Char *name) {
    
    f.open(name);        
    f << setprecision(PREC) << seed << "\n";     
    f.close();
    
}

        
// Helper functions from NR3, p. 635.
inline Complex* Cmplx(VecDoub &d) {return (Complex *)&d[0];}
inline Complex* Cmplx(Doub *d) {return (Complex *)d;}

void Gt_update(Int &k, MatDoub &data, 
    VecDoub &Gt00, VecDoub &Gt01, VecDoub &Gt11, VecDoub &Gt02, VecDoub &Gt12, 
    VecDoub &Gt22, VecDoub &Gt03) {
            
    Gt00[k] = abs(Cmplx(data[0])[0]) * abs(Cmplx(data[0])[0]); 
    Gt01[k] = abs(Cmplx(data[0])[1]) * abs(Cmplx(data[0])[1]); 
    Gt11[k] = abs(Cmplx(data[1])[1]) * abs(Cmplx(data[1])[1]); 
    Gt02[k] = abs(Cmplx(data[0])[2]) * abs(Cmplx(data[0])[2]); 
    Gt12[k] = abs(Cmplx(data[1])[2]) * abs(Cmplx(data[1])[2]); 
    Gt22[k] = abs(Cmplx(data[2])[2]) * abs(Cmplx(data[2])[2]); 
    Gt03[k] = abs(Cmplx(data[0])[3]) * abs(Cmplx(data[0])[3]);

}


// Construct the adjacency matrix as a sparse matrix.
NRsparseMat adjacency(Int &z, Int &L, Int &N, Int &nvals) { 
    
    Int i, j, k, l, site, pad = z;

    NRsparseMat S = NRsparseMat(N, N, nvals); 
    l = 0;  
    VecInt n(pad, 0);
    
    if (z == 2) {

        for (i = 0; i < N; ++i) {
              
            site = i;                 

            n[0] = (N+i-1) % N;
            n[1] = (N+1) % N;

            sort(n);

            for (k = 0; k < pad; ++k) {
                
	       S.val[l] = 1;
	       S.row_ind[l] = n[k];
	       ++l;
	           
            }

            S.col_ptr[site+1] = l;
                
        }   

        return S; 
               
    }
    
    if (z == 4) {
            
        for (i = 0; i < L; ++i) {
            for (j = 0; j < L; ++j) { 
              
                site = i * L + j;                 

                n[0] = i * L + ((L+j-1) % L);
                n[1] = i * L + ((j+1) % L);
                n[2] = ((L+i-1) % L) * L + j;
                n[3] = ((i+1) % L) * L + j;

                sort(n);

                for (k = 0; k < pad; ++k) {
                    
	           S.val[l] = 1;
	           S.row_ind[l] = n[k];
	           ++l;
	           
                }

                S.col_ptr[site+1] = l;
                
            }
        }   

        return S;
    
    }
}


// Construct the dispersal matrix as a sparse matrix.
NRsparseMat dispersal(Doub &coupling, Int &z, Int &L, Int &N, Int &nvals) { 
    
    Int i, j, k, l, site;
    Int pad = z;
  
    NRsparseMat S = NRsparseMat(N, N, nvals); 
    l = 0;
    VecInt n(pad+1, 0);

    if (z == 2) {
        
        for (i = 0; i < N; ++i) {
            
            site = i;
                    
            n[0] = (N+i-1) % N;
            n[1] = (i+1) % N;

            n[2] = site;

            sort(n);

            for (k = 0; k < pad+1; ++k) {
        	   S.val[l] = (n[k] == site) ? 1 - z * coupling : coupling;
        	   S.row_ind[l] = n[k];
        	   ++l;
            }
    
            S.col_ptr[site+1] = l;
                
        }

        return S;
        
    }
        
    if (z == 4) {
        
        for (i = 0; i < L; ++i) {
            for (j = 0; j < L; ++j) {
            
                site = i * L + j;
                    
                n[0] = i * L + ((L+j-1) % L);
                n[1] = i * L + ((j+1) % L);
                n[2] = ((L+i-1) % L) * L + j;
                n[3] = ((i+1) % L) * L + j;

                n[4] = site;

                sort(n);

                for (k = 0; k < pad+1; ++k) {
        	       S.val[l] = (n[k] == site) ? 1 - z * coupling : coupling;
        	       S.row_ind[l] = n[k];
        	       ++l;
                }
    
                S.col_ptr[site+1] = l;
                
            }
        }

        return S;
        
    }
  
}


// Initialize the population configuration.
VecDoub x_init(Doub &mean, Doub &std, Int &N, SNormaldev &SNormaldev) {
    
    VecDoub x(N, 0.0);
      
    for (Int i = 0; i < N; ++i) x[i] = mean + std * SNormaldev.dev();
    
    return x;
    
}


// Right multiplication of a sparse matrix by a vector.
void ax(NRsparseMat &S, VecDoub &x) {

    Int j, i;
    VecDoub y(S.ncols, 0.0);
            
    for (j = 0; j < S.ncols; j++) {
        
        for (i = S.col_ptr[j]; i < S.col_ptr[j+1]; i++) 
    
            y[S.row_ind[i]] += S.val[i] * x[j];
      
    }
    
    for (j = 0; j < S.ncols; j++) x[j] = y[j];
  
}


void ax(NRsparseMat &S, VecDoub &x, VecDoub &y) {

    Int j, i;
    y.assign(S.ncols, 0.0);
            
    for (j = 0; j < S.ncols; j++) {
        
        for (i = S.col_ptr[j]; i < S.col_ptr[j+1]; i++) 
    
            y[S.row_ind[i]] += S.val[i] * x[j];
      
    }
  
}


// Mean value of a vector.
Doub m(VecDoub &x, Int &N) {
  
    Doub m = 0.0;

    for (Int i = 0; i < N; ++i) m += x[i];

    return m / N;
    
}


// Energy.
Doub e(VecDoub &x, VecDoub &y, NRsparseMat &A, Int &N) {
    
    Doub e = 0.0;

    ax(A, x, y);
  
    for (Int i = 0; i < N; ++i) e += -x[i] * y[i];
    
    return e / N;
  
}


// Apply the Ricker to each element of a vector "x".
void ricker(VecDoub &x, Doub &lam, Doub &R, Int &N, SNormaldev &SNormaldev) {
  
    Doub mean;
    
    for (Int i = 0; i < N; ++i) {
        
        mean = x[i];
        
        if (mean > TOL) {
            
            mean = mean * exp(R * (1.0 - mean));
            mean = (mean > TOL) ? mean : 0.0;
            mean += SNormaldev.dev() * lam * mean;
            x[i] = (mean > TOL) ? mean : 0.0;
        
        } else {
            
            x[i] = 0.0;
        
        }
    }
}


void ricker(VecDoub &x, VecDoub &y, Doub &lam, Doub &R, Int &N, 
    SNormaldev &SNormaldev) {
  
    Doub mean;
    
    for (Int i = 0; i < N; ++i) {
        
        mean = x[i];
        
        if (mean > TOL) {
            
            mean = mean * exp(R * (1.0 - mean));
            mean = (mean > TOL) ? mean : 0.0;
            mean += SNormaldev.dev() * lam * mean;
            y[i] = (mean > TOL) ? mean : 0.0;
        
        } else {
            
            y[i] = 0.0;
        
        }
    }
}


/*------------------------------- Run --------------------------------*/

// Run the MC simulation for the DLR values designated for the current job.
Int mc(Int &d, Int &L, Int &lam_as_int, Char *suffix) {

    // Convert R_as_int to Doub.
    Doub lam = lam_as_int / Doub(LAM_MULT);
    
    cout << d << " " << L << " " << lam << endl;  
        
    // Derived parameters.
    Int z = 2 * d;
    Doub coupling = eps / z;    

    // Define the number of sites on the lattice.
    Int N = L * L;
    
    // Seed the random number generation.
    Uint seed = d + L + lam_as_int + clock();
    SNormaldev SNormaldev(seed);

    // Save the seed.
    sprintf(buf, "%sseed_%s", save_dir, suffix);
    save_to_file(seed, buf);
        
    // Generate the adjacency matrix.
    Int nval = z * N;
    NRsparseMat E = adjacency(z, L, N, nval);

    // Generate the disersal matrix.
    nval = (1 + z) * N;
    NRsparseMat M = dispersal(coupling, z, L, N, nval);

    // Initialize the population configuration.
    Doub std0 = sqrt(sig20);
    VecDoub x = x_init(mu0, std0, N, SNormaldev);
    
    // Iterate to burn in.
    cout << "starting burn-in at " 
         << ((Doub) clock() - start) / CLOCKS_PER_SEC << endl;
    
    Int i; 
    VecDoub y(N, 0.0);
    
    for (i = 0; i < num_sweeps_burnin; ++i) {

        // Two steps of the SCML equal one "sweep".
        ricker(x, lam, R, N, SNormaldev);
        ax(M, x, y);

        ricker(y, lam, R, N, SNormaldev);
        ax(M, y, x);
                
    }
    
    // Iterate asymptotics and gather statistics.
    cout << "starting asymptotics at " 
         << ((Doub) clock() - start) / CLOCKS_PER_SEC << endl;
    
    Doub xbar, ybar, xbar_mean, xbar2_mean, phi_mean, phi2_mean, value;
    VecDoub phi_list(num_stats, 0.0), e_list(num_stats, 0.0), Phi(N, 0.0), \
        speq(2 * L, 0.0), dummy(N, 0.0), \
        Gt00(num_stats, 0.0), Gt01(num_stats, 0.0), Gt11(num_stats, 0.0), \
        Gt02(num_stats, 0.0), Gt12(num_stats, 0.0), Gt22(num_stats, 0.0), \
        Gt03(num_stats, 0.0);
    MatDoub phi_array(num_x, N, 0.0);
    Int j, k = -1;       
                        
    for (i = 1; i < num_sweeps; ++i) {
        
        if (i >= x_record) 
            for (j = 0; j < N; ++j) phi_array[i - x_record][j] = (x[j] - y[j]) / 2;
        
        // Two steps of the SCML equal one "sweep".
        ricker(x, lam, R, N, SNormaldev);
        ax(M, x, y);

        ricker(y, dummy, lam, R, N, SNormaldev);
        ax(M, dummy, x);
         
        // Calculate running totals for xbar_mean and phi_mean.
        if (((i+1) % running_stats_record) == 0) { 
                                
            ybar = m(y, N);
            xbar = m(x, N);

            value = (xbar + ybar) / 2;
            xbar_mean += value;
            xbar2_mean += value * value;
            
            value = abs((xbar - ybar) / 2);
            phi_mean += value;
            phi2_mean += value * value;
            
            // Update lists of "phi", and "e".
            if (((i+1) % list_stats_record) == 0) { 
                                        
                phi_list[++k] = (xbar - ybar) / 2;
                
                for (j = 0; j < N; ++j) Phi[j] = (x[j] - y[j]) / 2;            
                //if (phi_list[k] < 0) for (j = 0; j < N; ++j) Phi[j] *= -1.0; 
                
                e_list[k] = e(Phi, dummy, E, N);
                
                // Only z == 4 case implemented here.   
                for (j = 0; j < N; ++j) dummy[j] = Phi[j];                 
                Doub *ptr = &dummy[0];
                MatDoub data(L, L, ptr);
    
                // The FFT of "data" rewrites "data" and "speq".  
                // (The initial values of "speq" do not affect the result.)
                rlft3(data, speq, 1);
    
                Gt_update(k, data, Gt00, Gt01, Gt11, Gt02, Gt12, Gt22, Gt03);                
            
            }

        } 
                    
    }

    // Save asymptotics.
    VecDoub running_stats(4, 0.0);
    running_stats[0] = xbar_mean / num_running_updates;
    running_stats[1] = xbar2_mean / num_running_updates;
    running_stats[2] = phi_mean / num_running_updates;
    running_stats[3] = phi2_mean / num_running_updates;
    
    sprintf(buf, "%srunning_stats_%s", save_dir, suffix);
    save_to_file(running_stats, buf);
    
    sprintf(buf, "%sphi_%s", save_dir, suffix);
    save_to_file(phi_list, buf);
    
    sprintf(buf, "%se_%s", save_dir, suffix);
    save_to_file(e_list, buf);

    sprintf(buf, "%sGt00_%s", save_dir, suffix);
    save_to_file(Gt00, buf);
    
    sprintf(buf, "%sGt01_%s", save_dir, suffix);
    save_to_file(Gt01, buf);
    
    sprintf(buf, "%sGt11_%s", save_dir, suffix);
    save_to_file(Gt11, buf);
    
    sprintf(buf, "%sGt02_%s", save_dir, suffix);
    save_to_file(Gt02, buf);
    
    sprintf(buf, "%sGt12_%s", save_dir, suffix);
    save_to_file(Gt12, buf);
    
    sprintf(buf, "%sGt22_%s", save_dir, suffix);
    save_to_file(Gt22, buf);
    
    sprintf(buf, "%sGt03_%s", save_dir, suffix);
    save_to_file(Gt03, buf);
    
    sprintf(buf, "%sx_%s", save_dir, suffix);
    save_to_file(phi_array, buf);
                      
    return 1;
    
}


int main(int argc, char** argv) {
         
    Int i, j, k;
    
    VecInt D_list(D_num, 0);    
    VecInt L_list(L_num, 0);

    for (i = 0; i < D_num; ++i) D_list[i] = D_low + i * D_inc; 
    for (i = 0; i < L_num; ++i) L_list[i] = L_low * int(pow(L_pow_inc, i));

    // Setup the lam list.
    VecDoub lam_list(lam_num, 0.0);
    Doub lam_delta = (lam_high - lam_low) / (lam_num - 1);
    for (i = 0; i < lam_num; ++i) lam_list[i] = lam_low + i * lam_delta;   
            
    // Get the command-line integer corresponding to the DLR values.
    Int item = atoi(argv[1]);
    cout << item << endl;
    
    // Only record DLR values for the first sarray job.
    if (item == 0) {
                
        char buf[CHAR_SIZE] = {0};
        sprintf(buf, "%sDLR.txt", save_dir);                
        f.open(buf);
        
    }

    // Generate the matrix of DLR values.
    MatInt DLR_list(D_num * L_num * lam_num + 1, 3, 0.0);
    Int index;
    
    for (i = 0; i < D_num; ++i) {
        for (j = 0; j < L_num; ++j) {
            for (k = 0; k < lam_num; ++k) {
            
                index = i * L_num * lam_num + j * lam_num + k;
                
                DLR_list[index][0] = D_list[i];
                DLR_list[index][1] = L_list[j];
                DLR_list[index][2] = int(lam_list[k] * LAM_MULT);
                
                // Only record DLR values for the first sarray job.
                if (item == 0) f << D_list[i] << " " 
                                 << L_list[j] << " " 
                                 << fixed << int(lam_list[k] * LAM_MULT) << endl;
            
   	    }
        }
    }

    // Only record DLR values for the first sarray job.
    if (item == 0) f.close();
                
    char suffix[CHAR_SIZE];    
    sprintf(suffix, "%dd_%d_%d.txt", \
        DLR_list[item][0], DLR_list[item][1], DLR_list[item][2]);

    Int return_value = \
        mc(DLR_list[item][0], DLR_list[item][1], DLR_list[item][2], suffix);

    cout << ((Doub) clock() - start) / CLOCKS_PER_SEC << endl;

    (return_value == 1) ? cout << "success" << endl : cout << "error" << endl;
    
    return 0;
    
}