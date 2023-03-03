target_sources(CUDATracer
	PRIVATE 
		${CMAKE_CURRENT_LIST_DIR}/BitmapWriter.cpp
		${CMAKE_CURRENT_LIST_DIR}/Camera.cpp
		${CMAKE_CURRENT_LIST_DIR}/Color.cpp
		${CMAKE_CURRENT_LIST_DIR}/Colorfunctions.cpp
		${CMAKE_CURRENT_LIST_DIR}/Direction.cpp
		${CMAKE_CURRENT_LIST_DIR}/Dualray.cpp
		${CMAKE_CURRENT_LIST_DIR}/Face.cpp
		${CMAKE_CURRENT_LIST_DIR}/Flatpoint.cpp
		${CMAKE_CURRENT_LIST_DIR}/GA_CPU_Tracer_glut.cpp
		${CMAKE_CURRENT_LIST_DIR}/GA_GPU_Tracer_glut.cpp
		${CMAKE_CURRENT_LIST_DIR}/LA_CPU_Tracer_glut.cpp
		${CMAKE_CURRENT_LIST_DIR}/Light.cpp
		${CMAKE_CURRENT_LIST_DIR}/Log.cpp
		${CMAKE_CURRENT_LIST_DIR}/main_glut.cpp
		${CMAKE_CURRENT_LIST_DIR}/Material.cpp
		${CMAKE_CURRENT_LIST_DIR}/MaterialBuilder.cpp
		${CMAKE_CURRENT_LIST_DIR}/Mesh.cpp
		${CMAKE_CURRENT_LIST_DIR}/Mesh2.cpp
		${CMAKE_CURRENT_LIST_DIR}/Plane.cpp
		${CMAKE_CURRENT_LIST_DIR}/Point.cpp
		${CMAKE_CURRENT_LIST_DIR}/Pointpair.cpp
		${CMAKE_CURRENT_LIST_DIR}/Ray.cpp
		${CMAKE_CURRENT_LIST_DIR}/Ray3D.cpp
		${CMAKE_CURRENT_LIST_DIR}/Sphere.cpp
		${CMAKE_CURRENT_LIST_DIR}/SplitScene.cpp
		${CMAKE_CURRENT_LIST_DIR}/timing.cpp
		${CMAKE_CURRENT_LIST_DIR}/Transformation.cpp
		${CMAKE_CURRENT_LIST_DIR}/Triangle.cpp
		${CMAKE_CURRENT_LIST_DIR}/VirtualTracer_glut.cpp
		${CMAKE_CURRENT_LIST_DIR}/kernel.cu
)

get_filename_component(SPPTracer ../ ABSOLUTE)

#target_include_directories(CUDATracer PUBLIC BoundSphere Meshreader)

