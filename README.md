# DUNE TMS
This is the project for studying the Temporary Muon Spectrometer as a part of the DUNE Near Detector system. 

It uses `edep-sim` output, which are stored at:

* Fourth production: `/pnfs/dune/persistent/ndmuonspect/TMSTarget/edep/`, using TMS as the only target, with a fixed TMS geometry (front 40 planes now in correct y), and no plane tilts. Used for pileup, flux studies, and so on

* Third production: `/pnfs/dune/persistent/users/marshalc/LArTMSProductionJun23withLArCV/edep/FHC/00m/00/`, using LAr as active target, with a slightly updated TMS geometry. Used for track matching studies by Faiza/Jeremy/Clarence

* Second production: `/pnfs/dune/scratch/users/marshalc/geomValHallLArTMS2/edep/0m/00/`, using LAr+TMS+Hall as active targets, with an updated more realistic cavern and cryostat, and updated TMS geometry.

* First alpha production: `/pnfs/dune/persistent/ndmuonspect/EDepSim_Sim`. Only LAr as active target. Some bugs in geometry. Only use these to reproduce Preliminary Design Report studies!

For more informatinon on running a production, see the [DUNE ND Production](https://github.com/DUNE/ND_Production) repo.

# Building the code
The repo is designed for minimal external dependencies. We use `ROOT` for some convenience classes, and input/output. There's dependency on `edep-sim` for the event and hit structure, and `toml` for config files.

## Setup and dependencies
The framework depends on `edep-sim`, `ROOT`, `CLHEP`, and `toml`. An example setup using mostly `ups` products at FNAL is provided in `setup.sh`: once the `edep-sim` `ups` product is set up, `ROOT` and `CLHEP` will also be appropriately set up.

## Building
### With Make
Once you have set your environment up, run `make`, which will get toml11, make the `src` directory and build the shared object (library) to `lib`, move onto the applications in `app` and build them into `bin`.

### With CMake (currentlty experimental)
To build with CMake, create a directory for the build (e.g. `mkdir -p build; cd build` inside dune-tms/), run `cmake ../`, and if all goes well you can them run `make`. Standard CMake options available, e.g. `-DCMAKE_INSTALL_PREFIX=/path/to/install/dir`, `-DCMAKE_CXX_FLAGS="-fpermissive -O999 -funroll-all-loops"`

## toml submodule dependency
We read parameter settings using `toml`. I particularly like [ToruNiina](https://github.com/ToruNiina/toml11/)'s repo, so have included it as a submodule here. The submodule setup is included in the default `make` recipe, so you shouldn't have to do anything else. If you for some reason want to run it on your own, do `git submodule init && git submodule update` to get the submodule.

The repo has a good readme, with fast and simple implementation, feel free to check it out and implement accordingly in the `TMS_Manager` class.

# Running the code
So far all the reconstruction happens through the `app/ConvertToTMSTree.cpp` driver application. In a nutshell, it takes `edep-sim` events, converts it into TMS formats, queries the geometry, runs reconstruction, and so on. There is a `TMS_TreeWriter` class which drives the `TTree` writing to file, and a `TMS_EventViewer` class which shows simple 2D event displays. The output `ROOT` file should contain (at least) two branches: one with reconstructed information and one with truth information. As of right October 2021, no hit level information is saved.

In the `app` directory you can also find some test executables for the reconstruction.

## Configuration options
All run-time configuration happens through `config/TMS_Default_Config.toml`. At time of writing, you can control reconstruction parameters (e.g. how many points do you need to form a "cluster" in DBSCAN, or what is the track reconstruction method), and application parameters (e.g. do we save detailed truth information, do we draw PDFs of hits) here. Have a look around and feel free to play!

# Directory structure
* `app` contains the example executables, linking to the TMS library
* `src` contains the TMS source files, like the track finder, event classes, true particle classes, and so on
* `scripts` contains simple scripts to run TMS studies in truth without reconstruction. Will produce output like `/pnfs/dune/persistent/ndmuonspect/FlatTrees`
* `utils` contains helper files mostly used for generating events with `edep-sim`. Most of the time you won't need these and they're mostly for documentation. You can also find these at `/pnfs/dune/persistent/ndmuonspect/Geometries`.

# Contact
* Liam O'Sullivan, [liam.osullivan@uni-mainz.de](mailto:liam.osullivan@uni-mainz.de)
* Gavin Davies, [gsdavies@phy.olemiss.edu](mailto:gsdavies@phy.olemiss.edu)
* Chris Marhsall, [chris.marshall@rochester.edu](mailto:chris.marshall@rochester.edu)
* Mathew Muether, [mathew.muether@wichita.edu](mailto:mathew.muether@wichita.edu)
* Clarence Wret (Inactive)
  
#nd\_muon\_spectrometer and #nd\_muon\_spectrometer\_code on [DUNE slack](https://dunescience.slack.com/)
