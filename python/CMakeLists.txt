cmake_minimum_required(VERSION 3.16)

project(PythonDCA)

file(GLOB SOURCES "bindings.cpp")

pybind11_add_module(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC DCA)

target_include_directories(${PROJECT_NAME} PUBLIC DCA)

IF(APPLE)
    SET_TARGET_PROPERTIES(${PROJECT_NAME}
        PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
        SUFFIX     ".so"
        )
ENDIF()

IF(UNIX AND NOT APPLE)
    set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".so")
endif()

IF(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".pyd")
endif()


set_target_properties(${PROJECT_NAME} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PythonDCA)

CONFIGURE_FILE(python_files/setup.py
   ${CMAKE_CURRENT_BINARY_DIR}/PythonDCA/setup.py
)
   
# Ditto for the setup file.
CONFIGURE_FILE(python_files/__init__.py
    ${CMAKE_CURRENT_BINARY_DIR}/PythonDCA/__init__.py
)
