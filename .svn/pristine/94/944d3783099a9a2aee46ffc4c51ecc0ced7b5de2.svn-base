INCLUDE_DIRECTORIES(.)

FILE(GLOB hdrs "*.h")
FILE(GLOB srcs_base "*.c*")

#evaluators
SET(eval_base evaluators/countforeg.cpp  evaluators/fogvh.cpp  evaluators/metaevaluator.cpp  evaluators/optflow.cpp  evaluators/optflowpm.cpp  evaluators/rgcolorhistogram.cpp  evaluators/surf.cpp  evaluators/xorcolor.cpp  evaluators/xor.cpp  evaluators/xorgrad.cpp evaluators/ramanan.cpp)
IF(gustereo3_FOUND)

SET(eval_st evaluators/stev.cpp  evaluators/stxorcolor.cpp  evaluators/stxor.cpp )
ENDIF()

IF(CUDA_FOUND)
SET(eval_gpu gpu/gpu_xor.cpp )
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_20,code=sm_20;-Xcompiler;-fPIC;)
cuda_compile(gpu_internal_OBJ gpu/gpu_xor.cu)
ENDIF()

#algorithms
FILE(GLOB srcs_algo ${srcs} "algorithms/ann/*.c*" "algorithms/cmaes/*.c*" "algorithms/condensation/*.c*" "algorithms/de/*.c*" "algorithms/macma/*.c*"  "algorithms/pso/*.c*")
FILE(GLOB hdrs_algo "algorithms/*.h" "algorithms/ann/*.h" "algorithms/cmaes/*.h" "algorithms/condensation/*.h" "algorithms/de/*.h" "algorithms/macma/*.h"   "algorithms/pso/*.h")

FILE( GLOB srcs ${srcs_base} ${srcs_eval} ${srcs_algo} ${eval_gpu} ${eval_st})
         

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} )

SET(hdrs ${hdrs} ${eval_base} )

ADD_LIBRARY(${PROJECT_NAME} ${srcs} ${hdrs} ${gpu_internal_OBJ})

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES # create *nix style library versions + symbolic links
	DEFINE_SYMBOL DSO_EXPORTS
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_SOVERSION}
	CLEAN_DIRECT_OUTPUT 1 										# allow creating static and shared libs without conflicts
	OUTPUT_NAME "${PROJECT_NAME}${PROJECT_DLLVERSION}" 	# avoid conflicts between library and binary target names
)

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${REQUIRED_LIBRARIES} )

INSTALL(TARGETS ${PROJECT_NAME}
	RUNTIME DESTINATION bin COMPONENT main
	LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE COMPONENT main
	ARCHIVE DESTINATION lib COMPONENT main)
MESSAGE(STATUS "hds= ${hdrs}")
INSTALL(FILES ${hdrs}
	DESTINATION include/${PROJECT_NAME}
	COMPONENT main)

INSTALL(FILES ${hdrs_evaluators}
	DESTINATION include/${PROJECT_NAME}/evaluators 
	COMPONENT main)

INSTALL(FILES ${hdrs_algo}
	DESTINATION include/${PROJECT_NAME}/algorithms 
	COMPONENT main)
