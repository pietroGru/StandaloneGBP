# StandaloneGBP absorber branch
Standalone Geant4 simulation for the Gamma Beam Profiler (GBP) detector

The application runs over geant4_10_07_p02.

To compile this application: mkdir a build directory, cd to it, run cmake (or ccmake) and build the project. The c++17 standard for gcc will be used. Then compile the project with
make install -j 16
where 16 is the number of parallel threads you want to use during compilation.
The make process will create an install dir where the binary and macros are present.

If you open the StandaloneGBP exec., it will automatically load the vis.mac macro for visualization. This macro takes care of opening the OpenGL window and initialize the run (geometry, phys. list, etc.).

In the vis mode, that is with the window opened, you can run gammaconversion.mac as an example of a macro where gaussian beam is used. This will generate 1000 photons and you will see the tracks on the screen.

If you want to perform any other macro, just drag & drop the macro over the application and let the simulation to finish. At the end, you will be prompted to press the 'Enter' key in order to close the window. As soon as you press Enter, the program will quit.
Alternatively, you can open a terminal window and pass the macro path to the exec. For instance:
./StandaloneGBP gaussianBeam.mac
