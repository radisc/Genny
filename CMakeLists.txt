cmake_minimum_required(VERSION 2.8.3)
project(Genny)

#FIND_PACKAGE( CPLEX REQUIRED)
find_package( CUDA )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -Xptxas -v,-gencode arch=compute_50,code=sm_50)


find_package (Threads)

include_directories( include	#${CPLEX_INCLUDE_DIRS}
				#${CONCERT_INCLUDE_DIR}
				src/ufl
				src/genetics	
				src/cuda)

#Executable
add_executable( Genny 		src/main.cpp
				src/ufl/ufl.h
				src/ufl/ufl.cpp
				src/genetics/genetics.h
				src/genetics/genetics.cpp)

cuda_add_library( cudaRO	src/cuda/cuda_functions.h
				src/cuda/cuda_functions.cu)

target_link_libraries( cudaRO 	${CUDA_LIBRARIES}
				${ILOCPLEX_LIBRARY}) 	
				#${CPLEX_LIBRARIES}
				#${CONCERT_LIBRARY})			

#cuda_add_library( cudaRO SHARED	src/cuda/cuda_functions.h
#					src/cuda/cuda_functions.cu)

#target_link_libraries( cudaRO	${CUDA_LIBRARIES})


#set_property(	TARGET Genny
# 		PROPERTY COMPILE_DEFINITIONS IL_STD
#  		)

target_link_libraries( Genny	cudaRO)
				#${ILOCPLEX_LIBRARY} 	
				#${CPLEX_LIBRARIES}
				#${CONCERT_LIBRARY})
