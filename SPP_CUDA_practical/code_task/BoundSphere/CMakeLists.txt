
add_library(BoundSphereLib 
	BoundSphereCreator.cpp
	Miniball.cpp 
	Miniball_b.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(BoundSphereLib PUBLIC OpenMP::OpenMP_CXX)
endif()


target_include_directories(BoundSphereLib PUBLIC ../CUDATracing)
