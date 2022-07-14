# Hybrid Plasma Simulation

This code is a space plasma simulation used for my research in the summer 2022 CSPAR/MSFC National Science Foundation Research Experience for Undergraduates. I build upon the code used by [Florinski et. al \(2016\)](https://doi.org/10.3847/0004-637X/826/2/197). 

## Simulation Methods

The code is a typical hybrid model, treating ions as individual particles and electrons as a massless fluid which maintains quasineutrality. 1 dimension is modeled in space while fields and velocities are modeled in 3 dimensions. It is a semi-implicit simulation, which uses an implicit method to advance the magnetic field and an explicit method to advance the electric field. 

## Acknowledgements

This work was supported by the NSF award AGS-1950831. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
