# Visualize a critical state of sync

By [Andrew Noble](http://two.ucdavis.edu/~andrewnoble)

## About

This repo contains the C++ and Python files used to simulate and animate a dynamical 2D Ising critical state in the collective synchronization of noisy coupled nonlinear two-cycle oscillators.  The end result is the animation at the top of  [this webpage](http://two.ucdavis.edu/~andrewnoble/research.html).

## Requirements

* g++
* Python (scipy, matplotlib, pylab)

## Usage

Clone the repo.
```
git clone https://github.com/andrewenoble/critical-sync.git
```
Run the Monte Carlo simulation written in C++ (using header files from Numerical Recipes 3rd ed).  200 1MB output files will be written to ```critical-sync/simulation_output```.  This may take a few minutes.  The place holder file ```simulation_output/m_0.txt``` will be overwritten.  
```
cd critical-sync/simulation
g++ asymp.cpp -O3 
./a.out
```
Generate the animation.  The existing animation ```critical_sync.mp4``` will be overwritten.
```
cd ../animation
python critical_sync_anim.py
```

## Acknowledgements

This work is support by an <a href="http://www.nsf.gov/awardsearch/showAward?AWD_ID=1344187&amp;HistoricalAwards=false">NSF
INSPIRE award</a> from the National Science Foundation.  
Details are published in <a href="http://www.nature.com/ncomms/2015/150408/ncomms7664/full/ncomms7664.html">Nature Communications</a>.
