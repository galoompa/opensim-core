


# Variable definitions
set(EXAMPLE_TARGET exampleVolumetricContactPlugin)
set(EXAMPLE_DIR ${OpenSim_SOURCE_DIR}/OpenSim/Examples/Plugins/VolumetricContactExample)
            
file(GLOB SOURCE_FILES ${EXAMPLE_DIR}/*.h ${EXAMPLE_DIR}/*.cpp)
add_library(${EXAMPLE_TARGET} SHARED ${SOURCE_FILES})

target_link_libraries(${EXAMPLE_TARGET} osimTools)
set_target_properties(${EXAMPLE_TARGET}
    PROPERTIES
    DEFINE_SYMBOL OSIMPLUGIN_EXPORTS
)
