cmake_minimum_required(VERSION 3.10.2) #change here (optional)
project(FastChemistryModel ) #change here

if (DEFINED ENV{WM_PROJECT})
    message("Using $ENV{WM_PROJECT}-$ENV{WM_PROJECT_VERSION}")
    set(WM_PATH, ${CMAKE_SOURCE_DIR})
else()
    message(FATAL_ERROR "OpenFOAM environment not set. Aborting.")
endif ()
message((${PROJECT_NAME}))

set(OpenFOAM_VERSION $ENV{WM_PROJECT_VERSION})
set(OpenFOAM_DIR $ENV{WM_PROJECT_DIR})
set(OpenFOAM_LIB_DIR $ENV{FOAM_LIBBIN})
set(OpenFOAM_SRC $ENV{FOAM_SRC})
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)  
set(CMAKE_BUILD_TYPE release)

#set(MKL_INCLUDE_DIR /opt/intel/oneapi/mkl/latest/include)
#set(MKL_LIB_DIR /opt/intel/oneapi/mkl/latest/lib/intel64)
#-g -fno-omit-frame-pointer
set(DEFINITIONS_COMPILE "-std=c++14 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP 
-DWM_LABEL_SIZE=32 -g -fno-omit-frame-pointer

-Wold-style-cast -Wnon-virtual-dtor 
-Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes   -DNoRepository 
-Wshadow -Wpedantic -Wsign-compare  -Wstrict-aliasing -Wnull-dereference -Wsequence-point
-Wreturn-local-addr -Wlogical-op  -Wduplicated-cond -Wduplicated-branches
-Wconversion
-Wall -Wextra -Wshadow -Winit-self -Wuninitialized -Wmaybe-uninitialized
-ftemplate-depth-100 -fPIC -O0 -mavx2  -mfma ")
#-fopt-info-vec-missed=vec_missed.log
add_definitions("${DEFINITIONS_COMPILE}")

include_directories(SYSTEM
        ${OpenFOAM_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude
        ${OpenFOAM_SRC}/MomentumTransportModels/compressible/lnInclude
        ${OpenFOAM_SRC}/ThermophysicalTransportModels/lnInclude
        ${OpenFOAM_SRC}/ThermophysicalTransportModels/fluidReactionThermo/lnInclude 
        ${OpenFOAM_SRC}/physicalProperties/lnInclude
        ${OpenFOAM_SRC}/thermophysicalModels/specie/lnInclude
        ${OpenFOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude 
        ${OpenFOAM_SRC}/thermophysicalModels/basic/lnInclude 
        ${OpenFOAM_SRC}/thermophysicalModels/chemistryModel/lnInclude
        ${OpenFOAM_SRC}/ODE/lnInclude
        ${OpenFOAM_SRC}/combustionModels/lnInclude
        ${OpenFOAM_SRC}/finiteVolume/lnInclude
        ${OpenFOAM_SRC}/finiteVolume/cfdTools
        ${OpenFOAM_SRC}/meshTools/lnInclude
        ${OpenFOAM_SRC}/sampling/lnInclude
        ${OpenFOAM_SRC}/OpenFOAM/lnInclude
        ${OpenFOAM_SRC}/OSspecific/POSIX/lnInclude

        )
include_directories( 

        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/asmCode       
        ${CMAKE_SOURCE_DIR}/ChemistryModel
        ${CMAKE_SOURCE_DIR}/ChemistryModel/Jacobian
        ${CMAKE_SOURCE_DIR}/dataBlock
        ${CMAKE_SOURCE_DIR}/dataBlock/simpleDataBlock
        ${CMAKE_SOURCE_DIR}/LUsolver
        ${CMAKE_SOURCE_DIR}/ODE/OptRodas34
        ${CMAKE_SOURCE_DIR}/ODE/OptRosenbrock34
        ${CMAKE_SOURCE_DIR}/ODE/OptSeulex
        ${CMAKE_SOURCE_DIR}/Reaction 
        ${CMAKE_SOURCE_DIR}/Reaction/updateReactionRate
        ${CMAKE_SOURCE_DIR}/Reaction/updateJacobianMatrix
        ${CMAKE_SOURCE_DIR}/Reaction/SRI
        ${CMAKE_SOURCE_DIR}/Reaction/Troe
        ${CMAKE_SOURCE_DIR}/Reaction/Lindemann
        ${CMAKE_SOURCE_DIR}/Thermo/idealGas 
        ${CMAKE_SOURCE_DIR}/Thermo
        ) #change here

link_directories(${OpenFOAM_LIB_DIR} ${OpenFOAM_LIB_DIR}/dummy ${OpenFOAM_LIB_DIR}/${PATH_LIB_OPENMPI} )
enable_language(ASM)
set_source_files_properties(
    asmCode/svml_d_gnuexp_data.S 
    asmCode/svml_d_gnuexp4_core_avx2.S 
    asmCode/svml_d_gnupow_data.S
    asmCode/svml_d_gnupow4_core_avx2.S
    asmCode/svml_d_gnulog_data.S
    asmCode/svml_d_gnulog4_core_avx2.S
    PROPERTIES 
    LANGUAGE ASM
)
add_library(
    ${PROJECT_NAME} SHARED 
    odeChemistrySolvers.C 
    asmCode/svml_d_gnuexp_data.S
    asmCode/svml_d_gnuexp4_core_avx2.S
    asmCode/svml_d_gnupow_data.S
    asmCode/svml_d_gnupow4_core_avx2.S
    asmCode/svml_d_gnulog_data.S
    asmCode/svml_d_gnulog4_core_avx2.S
    dataBlock/dataBlock.C 
    dataBlock/simpleDataBlock/simpleDataBlock.C
    Thermo/createIdealGasFromFoamDict.C
    Thermo/idealGas/idealGas.C
)

target_link_libraries(${PROJECT_NAME} 
fluidThermophysicalModels
specie
chemistryModel
ODE
combustionModels
reactionThermophysicalModels
momentumTransportModels
compressibleMomentumTransportModels
fluidReactionThermophysicalTransportModels
thermophysicalTransportModels
finiteVolume
fvModels
fvConstraints
meshTools
sampling
OpenFOAM dl Pstream
triSurface
)
