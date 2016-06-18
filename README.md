# Visualizing a critical state of sync

By [Andrew Noble](http://two.ucdavis.edu/~andrewnoble)

## About

This folder contains the C++ and Python files used to simulate and animate a dynamical 2D Ising critical state in the collective synchronization of noisy coupled nonlinear two-cycle oscillators.  The end result is the animation at the top of  [this webpage](http://two.ucdavis.edu/~andrewnoble/research.html).

## Requirements

* gcc
* Python 

## Get started

Clone the repo.
```
git clone https://github.com/andrewenoble/critical-sync.git
```
Run the Monte Carlo simulation written in C++.  Output files will be written to ```critical-sync/simulation_output```.  The place holder ```simulation_output/m_0.txt``` will be overwritten.
```
cd critical-sync/simulation
g++ asymp.cpp -O3 
```
Generate the animation.  The existing animation ```critical_sync.mp4``` will be overwritten.
```
cd ../animation
python critical_sync_anim
```

## Acknowledgements

These are preliminary results emerging from a collaboration with [Alan Hastings](http://two.ucdavis.edu/~me), [Jonathan Machta](http://people.umass.edu/machta), [Ottar Bjornstad](http://ento.psu.edu/directory/onb), and [Bryan Grenfell](https://www.princeton.edu/step/people/faculty/bryan-grenfell).
