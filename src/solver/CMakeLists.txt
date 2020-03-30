#Todo: add other solvers

if (NOT ENABLE_CUDA)
    if (SOLVER STREQUAL "hlle")
        set(SELECTED_SWE_SOLVER_HEADER HLLEFun.hpp)
    elseif (SOLVER STREQUAL "fwave")
        set(SELECTED_SWE_SOLVER_HEADER FWaveVec.hpp)
    elseif(SOLVER STREQUAL "rusanov")
    elseif(SOLVER STREQUAL "augRie")
        if(CUDA_ENABLED)
            set(SELECTED_SWE_SOLVER_HEADER AugRieCUDA.h)
        else()
            set(SELECTED_SWE_SOLVER_HEADER AugRieFun.hpp)
        endif()
    else()
        message (WARNING "${SOLVER} not recognized!")
    endif()
endif()

target_sources(swe 
    PUBLIC
    ${SELECTED_SWE_SOLVER_HEADER}
)
