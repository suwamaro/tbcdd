add_library(prog_base OBJECT
parameters_util.cpp
cpptoml_helper.cpp
)

add_executable(tbc tbc.cpp $<TARGET_OBJECTS:prog_base>)
target_link_libraries(tbc ${LIBS})
