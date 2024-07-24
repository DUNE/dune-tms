# Maintainer: Liam

# Declare the project to meson
project('dune-tms', 'cpp')

# Declare a subproject for a toml11 library, will be cloned and built on configure
toml = subproject('toml11')

# Declare the dependencies required to build, ROOT is a bit funny with meson
toml_dep = toml.get_variable('toml11_dep')
geant4_dep = dependency('Geant4', version : '>=4.10.6')
edep_dep = dependency('EDepSim')
root_dep = dependency('ROOT', version : '>=5.32.00')

#root_inc = include_directories('/uni-mainz.de/homes/losulliv/nongentoo-soft/root/Linux/include/')
#root_lib = shared_library('/uni-mainz.de/homes/losulliv/nongentoo-soft/root/Linux/lib/')
#root_dep = declare_dependency(link_with : root_lib, include_directories : root_inc)

deps = [root_dep, geant4_dep, edep_dep, toml_dep]

source_files = ['app/ConvertToTMSTree.cpp', 'src/BField_Handler.cpp', 'src/TMS_Bar.cpp',
                'src/TMS_Event.cpp', 'src/TMS_EventViewer.cpp', 'src/TMS_Hit.cpp', 'src/TMS_Kalman.cpp',
                'src/TMS_Manager.cpp', 'src/TMS_ReadoutTreeWriter.cpp', 'src/TMS_Readout_Manager.cpp',
                'src/TMS_Reco.cpp', 'src/TMS_Track.cpp', 'src/TMS_TreeWriter.cpp', 'src/TMS_TrueHit.cpp',
                'src/TMS_TrueParticle.cpp', 'src/TMS_Utils.cpp']
header_files = ['src/BField_Handler.h', 'src/BetheBloch.h', 'src/Material.h', 'src/MultipleScattering.h',
                'src/TMS_Bar.h', 'src/TMS_Constants.h', 'src/TMS_DBScan.h', 'src/TMS_Event.h',
                'src/TMS_EventViewer.h', 'src/TMS_Geom.h', 'src/TMS_Hit.h', 'src/TMS_Kalman.h',
                'src/TMS_Manager.h', 'src/TMS_ReadoutTreeWriter.h', 'src/TMS_Readout_Manager.h',
                'src/TMS_Reco.h', 'src/TMS_Track.h', 'src/TMS_TreeWriter.h', 'src/TMS_TrueHit.h',
                'src/TMS_TrueParticle.h', 'src/TMS_Utils.h']

executable('ConvertToTMSTree', sources : [source_files, header_files], dependencies : deps)

message('')
message('This meson build script is experimental. Report any issues to Liam <liam.osullivan@uni-mainz.de>')
message('')
