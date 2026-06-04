cmake_minimum_required(VERSION 3.10)
project(test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(test test.cpp)

target_link_libraries(test
	PRIVATE
		Integrator
		OpenMP::OpenMP_CXX
)

target_compile_options(test PRIVATE
	$<$<CONFIG:Release>:-O2>
	$<$<CONFIG:RelWithDebInfo>:-O2>
)
add_subdirectory(Integrator)