cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set options for this project.
set (PROJECT_NAME "pfsp" CXX)
project (${PROJECT_NAME})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# find glog
add_subdirectory(glog)

# brew install llvm libomp
# https://stackoverflow.com/questions/16182259/how-to-install-the-c-binding-for-zeromq-on-mac-os-x/21388687
if(APPLE)
    set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -stdlib=libc++")
endif()

#set(CMAKE_C_FLAGS "${CMAKE_CFLAGS} -fopenmp")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG")
if(MSVC)
    set(CMAKE_CXX_FLAGS "-DIL_STD /std:c++11 /permisive-")
else()
    set(CMAKE_CXX_FLAGS "-g -ggdb -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions -std=c++14")
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native  -O3")  # -fopenmp
set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_BUILD_TYPE Debug)

set (PROJECT_SOURCE_DECLARATION_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
set (PROJECT_SOURCE_DEFINITION_DIRECTORY ${PROJECT_SOURCE_DIR}/src/cpp)
set (MAIN_FILE ${PROJECT_SOURCE_DEFINITION_DIRECTORY}/main.cpp)

# Set options used by "CMake".
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

include_directories(src/cpp/controller/include)
include_directories(src/cpp/problem/branchandbound/include)
include_directories(src/cpp/problem/heuristic/include)
include_directories(src/cpp/problem/include)
include_directories(src/cpp/util/include)

file(GLOB APP_SOURCES */*.cpp)
foreach (testsourcefile ${APP_SOURCES})
    get_filename_component(testname ${testsourcefile} NAME_WLE)
    get_filename_component(dirname ${testsourcefile} DIRECTORY)
    file(GLOB dir_src CONFIGURE_DEPENDS
            "${dirname}/*.h"
            "${dirname}/*.cpp"
            )
    message("${testname}.cpp | ${dir_src}")
    add_executable("${testname}.cpp" ${dir_src})
endforeach (testsourcefile ${APP_SOURCES})

add_executable(${PROJECT_NAME}
        src/cpp/controller/PFSPFileController.cpp
        src/cpp/problem/deterministic/heuristic/include/HeuristicCommons.h
        src/cpp/problem/deterministic/heuristic/include/JohnsonSolver.h
        src/cpp/problem/deterministic/heuristic/include/NEHSolver.h
        src/cpp/problem/deterministic/heuristic/HeuristicCommons.cpp
        src/cpp/problem/deterministic/heuristic/JohnsonSolver.cpp
        src/cpp/problem/deterministic/heuristic/NEHSolver.cpp
        src/cpp/problem/deterministic/PFSInstance.h
        src/cpp/problem/deterministic/PFSProblem.h
        src/cpp/problem/PFSSolution.h
        src/cpp/problem/deterministic/PFSInstance.cpp
        src/cpp/problem/deterministic/PFSProblem.cpp
        src/cpp/problem/PFSSolution.cpp
        src/cpp/util/include/CollectionUtil.h
        src/cpp/util/include/EnumUtil.h
        src/cpp/util/include/ExecutionInfo.h
        src/cpp/util/include/FileUtil.h
        src/cpp/util/include/RandomUtil.h
        src/cpp/util/include/TimeDateUtil.h
        src/cpp/util/CollectionUtil.cpp
        src/cpp/util/EnumUtil.cpp
        src/cpp/util/ExecutionInfo.cpp
        src/cpp/util/FileUtil.cpp
        src/cpp/util/RandomUtil.cpp
        src/cpp/util/TimeDateUtil.cpp
        src/cpp/main.cpp
        src/cpp/util/NumericUtil.h
        src/cpp/controller/include/CommandLineInterfaceController.h
        src/cpp/controller/include/PFSPFileController.h
        src/cpp/controller/CommandLineInterfaceController.cpp
        src/cpp/util/NumericUtil.cpp
        src/cpp/util/include/SchedulingUtil.h
        src/cpp/problem/GlobalTypes.h
        src/cpp/problem/metaheuristic/common/Construction.cpp
        src/cpp/problem/metaheuristic/common/Construction.h
        src/cpp/problem/metaheuristic/common/LocalSearch.cpp
        src/cpp/problem/metaheuristic/common/LocalSearch.h
        src/cpp/problem/metaheuristic/common/Job.cpp
        src/cpp/problem/metaheuristic/common/Job.h
        src/cpp/problem/metaheuristic/common/SchedulingSolution.cpp
        src/cpp/problem/metaheuristic/common/SchedulingSolution.h
        src/cpp/problem/metaheuristic/common/RandNEHT.cpp
        src/cpp/problem/metaheuristic/common/RandNEHT.h
        src/cpp/problem/metaheuristic/common/Inputs.cpp
        src/cpp/problem/metaheuristic/common/Inputs.h
        src/cpp/problem/metaheuristic/common/Randomness.cpp
        src/cpp/problem/metaheuristic/common/Randomness.h
        src/cpp/problem/metaheuristic/grasp/GRASP.cpp
        src/cpp/problem/metaheuristic/grasp/GRASP.h
        src/cpp/problem/metaheuristic/common/ElapsedTime.cpp
        src/cpp/problem/metaheuristic/common/ElapsedTime.h
        src/cpp/problem/deterministic/branchandbound/DFSSearchStrategy.cpp
        src/cpp/problem/deterministic/branchandbound/DFSSearchStrategy.h
        src/cpp/problem/deterministic/branchandbound/SelectionStrategy.h
        src/cpp/util/random.h
        src/cpp/util/random.cpp
        src/cpp/problem/deterministic/branchandbound/PartialPermutation.h
        src/cpp/problem/robust/RobPFSInstance_Cmax.h
        src/cpp/util/singleton.h
        src/cpp/problem/robust/RobustBranchBound.h
        src/cpp/problem/robust/RobustBranchBound_Node.h
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Budget_Cmax_BranchBound.h
        src/cpp/problem/robust/PFSPScenario.h
        src/cpp/PFSP_Facade_Robust.h
        src/cpp/PFSP_Facade_Deterministic.h
        src/cpp/problem/PFSP_Parameters.h
        src/cpp/problem/metaheuristic/common/Test.cpp
        src/cpp/problem/metaheuristic/common/Test.h
        src/cpp/util/RandomStream.cpp src/cpp/util/RandomStream.h
        src/cpp/problem/deterministic/metaheuristic/GRASPSolver_PFSP_Cmax.cpp
        src/cpp/problem/deterministic/metaheuristic/GRASPSolver_PFSP_Cmax.h
        src/cpp/problem/metaheuristic/common/Outputs.cpp
        src/cpp/problem/metaheuristic/common/Outputs.h
        src/cpp/problem/robust/budget/PFSPBudgetScenario.h
        src/cpp/problem/robust/budget/PFSP_Cmax_Budget_WorstCase.h
        src/cpp/problem/robust/budget/metaheuristic/cmax/PFSP_Cmax_Budget_UpperBound_Alt.cpp
        src/cpp/problem/robust/budget/metaheuristic/cmax/PFSP_Cmax_Budget_UpperBound_Alt.h
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Cmax_Budget_LowerBound.cpp
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Cmax_Budget_LowerBound.h
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Cmax_Budget_NodeData.cpp
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Cmax_Budget_NodeData.h
        src/cpp/problem/robust/budget/branchandbound/cmax/PFSP_Budget_Cmax_BranchBound.h
        src/cpp/problem/robust/budget/metaheuristic/wct/GRASPSolver_RobPFSP_WCT.cpp
        src/cpp/problem/robust/budget/metaheuristic/wct/GRASPSolver_RobPFSP_WCT.h
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_LowerBound.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_LowerBound.h
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_NodeData.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_NodeData.h
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_Budget_WCT_BranchBound.h
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_Hybrid_BranchBound.h
        src/cpp/problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_Hybrid_BranchBound.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_DFSNodeHandler.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_DFSNodeHandler.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_BranchInfo.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_BranchInfo.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_Brancher.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_Brancher.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_SimpleNodeInfo.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/WilsonPFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/WilsonPFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TS2PFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TS2PFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TS3PFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TS3PFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/WST2PFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/WST2PFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TBAPFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/TBAPFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/LiaoYouPFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/LiaoYouPFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/MannePFSP_WCT_Model.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/MannePFSP_WCT_Model.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_ModelStrategy.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_ModelSelector.h
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_CPLEXEnvironment.cpp
        src/cpp/problem/robust/budget/branchandbound/wct/cplex/PFSP_WCT_CPLEXEnvironment.h
        src/cpp/problem/robust/budget/branchandbound/BB_OutputUtil.h
        src/cpp/problem/robust/budget/metaheuristic/cmax/GRASPSolver_RobPFSP_Cmax.cpp
        src/cpp/problem/robust/budget/metaheuristic/cmax/GRASPSolver_RobPFSP_Cmax.h
        src/cpp/problem/metaheuristic/common/OutputUtil.h
        src/cpp/util/precision.h src/cpp/util/precision.h
        src/cpp/problem/robust/RobPFSInstance.h
        src/cpp/problem/robust/RobPFSInstance_WCT.h
        src/cpp/problem/robust/budget/PFSP_WCT_Budget_WorstCase.cpp
        src/cpp/problem/robust/budget/PFSP_WCT_Budget_WorstCase.h 
        src/cpp/problem/robust/budget/metaheuristic/wct/GRASPSolver_RobPFSP_WCT.cpp 
        src/cpp/problem/robust/budget/metaheuristic/wct/GRASPSolver_RobPFSP_WCT.h 
        src/cpp/PFSP_Facade_Robust.cpp 
        src/cpp/PFSP_Facade_Deterministic.cpp)

#set(Boost_NO_SYSTEM_PATHS TRUE)
if (Boost_NO_SYSTEM_PATHS)
    set(BOOST_ROOT "/opt/boost")
    set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
    set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
endif (Boost_NO_SYSTEM_PATHS)

find_package(Boost 1.65.0 COMPONENTS system filesystem chrono timer program_options random date_time REQUIRED)  # thread

if (Boost_FOUND)
    # From the offical documentation:
    # Add include directories to the build. [...] If the SYSTEM option is given,
    # the compiler will be told the directories are meant as system include
    # directories on some platforms (signalling this setting might achieve effects
    # such as the compiler skipping warnings [...])."
    #include_directories (SYSTEM ${Boost_INCLUDE_DIR})
    #ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
    ADD_DEFINITIONS( "-DHAS_BOOST" )

    # From the offical documentation:
    # "Specify directories in which the linker will look for libraries. [...] Note
    # that this command is rarely necessary. Library locations returned by
    # find_package() and find_library() are absolute paths. Pass these absolute
    # library file paths directly to the target_link_libraries() command. CMake
    # will ensure the linker finds them."
    #link_directories (${Boost_LIBRARY_DIRS})
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

# using pkg-config to configure uuid
#set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(PkgConfig REQUIRED)
#find_package(Threads REQUIRED)

#pkg_search_module(UUID REQUIRED uuid)
# On MacOSX we have:
# $ nm -g /usr/lib/libSystem.dylib | grep uuid_generate
# 000b3aeb T _uuid_generate
# 0003e67e T _uuid_generate_random
# 000b37a1 T _uuid_generate_time
if(APPLE)
    set(UUID_LIBRARY_VAR System)
    include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
    link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
else()
    # Linux type:
    set(UUID_LIBRARY_VAR uuid)
endif()

find_library(UUID_LIBRARY
        NAMES ${UUID_LIBRARY_VAR}
        PATHS /lib /usr/lib /usr/local/lib
        )

find_package(OpenMP)
if(MSVC)
    set(OpenMP_LIBRARY_DIR "C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29910/lib/x64")
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Tell CMake where to find the OpenMP libraries.
    link_directories(${OpenMP_LIBRARY_DIR})
endif()

######## Cplex Dependency
find_package(CPLEX REQUIRED)
if (CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_LIBRARY})
endif ()

if (CPLEX_CONCERT_FOUND)
    include_directories(${CPLEX_CONCERT_INCLUDE_DIR})
    link_directories(${CPLEX_CONCERT_LIBRARY})
endif ()

if (CPLEX_ILOCPLEX_FOUND)
    include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_ILOCPLEX_LIBRARY})
endif ()

#include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ${Boost_INCLUDE_DIRS} )
#link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} )

set( CommonLibs
        m
        ilocplex
        cplex-library
        cplex-concert
        #pthread
        dl  # since CPLEX 12.8
	#boost_thread.so.1.65.1
	#boost_system.so.1.65.1
        )

# "Add an executable to the project using the specified source files."
#add_executable (${PROJECT_NAME} ${MAIN_FILE})
# "Link a target to given libraries."
if(APPLE)  # cplex-library cplex-concert ilocplex
    target_link_libraries (${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES} ${UUID_LIBRARIES} ${MATH_LIBRARY} ${CMAKE_DL_LIBS} glog::glog OpenMP::OpenMP_CXX ${CommonLibs})
elseif(MSVC)
    target_link_libraries (${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES} ${UUID_LIBRARIES} ${MATH_LIBRARY} ${CMAKE_DL_LIBS} ${OpenMP_CXX_FLAGS} glog::glog ${CommonLibs})
else()
    target_link_libraries (${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES} ${UUID_LIBRARIES} ${MATH_LIBRARY} ${CMAKE_DL_LIBS} glog::glog OpenMP::OpenMP_CXX rt ${CommonLibs})
endif()
