# DUNE TMS
This is the project for studying the Temporary Muon Spectrometer as a part of the DUNE Near Detector system. 

It uses `edep-sim` output, which are stored at:

* Second production: `/pnfs/dune/scratch/users/marshalc/geomValHallLArTMS2/edep/0m/00/`, using LAr+TMS+Hall as active target, with an updated more realistic cavern and cryostat, and updated TMS geometry.

* First alpha production: `/pnfs/dune/persistent/ndmuonspect/EDepSim_Sim`. Only LAr as active target. Some bugs in geometry. Only use these to reproduce Preliminary Design Report studies!

For more informatinon on running a production, see the https://github.com/DUNE/ND_Production repository.

# Running the code
So far all the interfacing happens through the `app/ConvertToTMSTree.cpp` driver application. In a nutshell, it takes `edep-sim` events, converts it into TMS formats, queries the geometry, runs reconstruction, and so on. This will hopefully be more elegant soon!

In the `app` directory you can also find some test executables for the reconstruction.

# Setup and dependencies
The framework depends on `edep-sim`, `ROOT`, and `CLHEP`. An example setup using mostly `ups` products is provided in `setup.sh`.

Once you have set your environment up, run `make`, which will make the `src` directory and build the shared object (library) to `lib`, move onto the applications in `app` and build them into `bin`.

## toml submodule dependency
We read parameter settings using `toml`. I particularly like the https://github.com/ToruNiina/toml11/ repo, so have included it as a submodule here. The repo has a good readme, with fast and simple implementation, feel free to check it out and implement accordingly in the `TMS\_Manager` class

# Directory structure
* `app` contains the example executables, linking to the TMS library
* `src` contains the TMS source files, like the track finder, event classes, true particle classes, and so on
* `scripts` contains simple scripts to run TMS studies in truth without reconstruction. Will produce output like `/pnfs/dune/persistent/ndmuonspect/FlatTrees`
* `utils` contains helper files mostly used for generating events with `edep-sim`. Most of the time you won't need these and they're mostly for documentation. You can also find these at `/pnfs/dune/persistent/ndmuonspect/Geometries`.

# Contact
* Clarence Wret, [c.wret@rochester.edu](mailto:c.wret@rochester.edu)
* Gavin Davies, [gsdavies@phy.olemiss.edu](mailto:gsdavies@phy.olemiss.edu)
* Chris Marhsall, [chris.marshall@rochester.edu](mailto:chris.marshall@rochester.edu)
* Mathew Muether, [mathew.muether@wichita.edu](mailto:mathew.muether@wichita.edu)

#nd\_muon\_spectrometer on dunescience.slack.com
