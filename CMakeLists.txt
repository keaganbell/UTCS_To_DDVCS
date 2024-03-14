#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8...3.19 FATAL_ERROR)
##project(AnaEx02)
project(tcs)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/shared/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/shared/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/shared/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
##add_executable(AnaEx02 AnaEx02.cc ${sources} ${headers})
##target_link_libraries(AnaEx02 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
add_executable(tcs tcs.cc ${sources} ${headers})
target_link_libraries(tcs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build AnaEx02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
##set(AnaEx02_SCRIPTS
##    AnaEx02.in run.mac init_vis.mac vis.mac
##  )
set(tcs_STUFF
  run.mac vis.mac beam_definition.txt field_param.txt random_seeds.dat
  stacking_control.txt field/g2p_hallbfield.dat tcs_gen.data tcs_gen.kin_data
  init_vis.mac nps_module.inp howto.txt
  tcs_gdmls/beam_pipe.gdml tcs_gdmls/cps.gdml tcs_gdmls/magnet.gdml
  tcs_gdmls/simple_magnet.gdml 
  tcs_gdmls/scattering_chamber_hallac_may4_2023.gdml tcs_gdmls/target_hallac.gdml	
  tcs_gdmls/target.gdml tcs_gdmls/calorimeter.gdml tcs_gdmls/lheshield.gdml
  tcs_gdmls/tcs_setup.gdml tcs_gdmls/carbon_target.gdml
  tcs_gdmls/ln2shield.gdml tcs_gdmls/scattering_chamber.gdml
  tcs_gdmls/pointer_referenced/tracker1_ref.gdml
  tcs_gdmls/pointer_referenced/tracker2_ref.gdml
  tcs_gdmls/pointer_referenced/tracker3_ref.gdml
  tcs_gdmls/pointer_referenced/tracker4_ref.gdml

   tcs_gdmls/scattering_chamber_hallac_may4_2023_test.gdml
)

##foreach(_script ${AnaEx02_STUFF})
foreach(_script ${tcs_STUFF})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
##install(TARGETS AnaEx02 DESTINATION bin)
install(TARGETS tcs DESTINATION bin)

