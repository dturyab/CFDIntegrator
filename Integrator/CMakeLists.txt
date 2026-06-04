add_library(Integrator STATIC)

file(GLOB INTEGRATOR_SRC "src/*.cpp")
file(GLOB INTEGRATOR_HS "include/*.h")

target_sources(Integrator
	PRIVATE
		${INTEGRATOR_SRC}
)

target_include_directories(Integrator
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Integrator
	PUBLIC
		OpenMP::OpenMP_CXX
)

target_compile_options(Integrator PRIVATE
	$<$<CONFIG:Release>:-O2>
    	$<$<CONFIG:RelWithDebInfo>:-O2>
)
