# AutoTomoAlign

We provide MATLAB sample Scripts and Functions for automatic 3D alignment and tomographic reconstruction of phase-contrast projection data.

The presented scripts support 2D and 3D parallel beam geometries with flexible definition of the sample/detector relative orientation.

The three-dimensional projection alignment is done in two different stages: First, a cross-correlation alignment in the projections gradient domain aims to correct relative translational movements between adjacent projections. The refinement of the translational parameters and the estimation of the angular projection errors is done by solving an optimisation problem via the Levenberg-Marquardt algorithm.

For phase-contrast projection data, coming for example from Coherent Diffraction Imaging (CDI) techniques, the proposed algorithm automatically corrects for linear background phase terms. 

Our modified SIRT reconstruction algorithm accounts for phase-wrapping in the projections under reasonable limitations (for real data applications).


## Documentation / samples

A detalied description of the developed algorithm is presented in (*ADD REFERENCE TO PUBLICATION*) where the results obtained from the sample script AutoTomoAlign/Scripts_1/Run_SheppLogan_Script.m are published.

Additional test Sripts are provided


## Installation instructions

The tomographic reconstruction algorithm implemented in this release makes use of the ASTRA toolbox. Please follow the instructions provided [here](https://github.com/astra-toolbox/astra-toolbox) depending on your operating system.

Please ensure that your current MATLAB installation includes the Imaging Processing Toolbox.

In MATLAB, simply run the test scripts in Scripts/ or modify the simulation parameters by editting the source script:

```edit Run_SheppLogan_Script```


## References

If you pretend to cite any part of the developed algorithms, we would apreciate it if you would refer to the following paper:
(*ADD REFERENCE TO PUBLICATION*).

## License

All files with exception of those in the External/ folder are open source under the [3-Clause BSD license](https://github.com/tiagoapcramos/AutoTomoAlign/blob/master/LICENSE.txt).

## Contacts

name:  Tiago Ramos <br />
email: tiagoj@dtu.dk


name:  Jakob Sauer Jørgensen<br />
email: jakj@dtu.dk


name:  Jens Wenzel Andreasen<br />
email: jewa@dtu.dk



