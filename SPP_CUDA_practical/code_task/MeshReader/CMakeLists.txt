
add_library(MeshReaderLib list.cpp
MeshFileReader.cpp
objLoader.cpp
obj_parser.cpp
rtmLoader.cpp
string_extra.cpp
string_util.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(BoundSphereLib PUBLIC OpenMP::OpenMP_CXX)
endif()


target_include_directories(MeshReaderLib PUBLIC ../CUDATracing)
