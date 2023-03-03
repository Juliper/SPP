cmake_minimum_required(VERSION 3.13)

# set the project name and version
project(CUDATracer VERSION 1.5)

# Initialize for C and C++ languages.
project(cmake_and_cuda)
# Also, enable CUDA language support.
enable_language(CUDA)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#list(APPEND EXTRA_LIBS BoundSphereLib MeshReaderLib)
#list(APPEND EXTRA_INCLUDES ${PROJECT_SOURCE_DIR}/BoundSphere 
#	${PROJECT_SOURCE_DIR}/MeshReader
#	${PROJECT_SOURCE_DIR}/CUDATracing)

# add the libraries
add_subdirectory(BoundSphere)
add_subdirectory(MeshReader)

# add the executable
add_executable(CUDATracer)
add_subdirectory(CUDATracing)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(CUDATracer PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(CUDATracer PUBLIC BoundSphereLib MeshReaderLib)

# add the binary tree to the search path for include files
target_include_directories(CUDATracer PRIVATE CUDATracing BoundSphere MeshReader )

# Get the include directories for the target.
get_target_property(CUDATINCLUDES CUDATracer INCLUDE_DIRECTORIES)
get_target_property(MESHRINLCUDES MeshReaderLib INCLUDE_DIRECTORIES)
get_target_property(BOUNDINCLUDES BoundSphereLib INCLUDE_DIRECTORIES)

message(STATUS "Cudatracer includes ${CUDATINCLUDES}")
message(STATUS "MeshReder includes ${MESHRINLCUDES}")
message(STATUS "BoundSphereLib includes ${BOUNDINCLUDES}")
