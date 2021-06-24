find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

if(SWIG_VERSION GREATER_EQUAL 4)
    set(CMAKE_SWIG_FLAGS "-doxygen")
endif()
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
include_directories(../../include)
#set_property(SOURCE foo.i PROPERTY SWIG_MODULE_NAME pytitanlibcustom)


# Add swig module
#set(CMAKE_SWIG_FLAGS -interface _pytitanlibcustom)
swig_add_library(titanlibcustom TYPE SHARED LANGUAGE python SOURCES ../${PROJECT_NAME}.i)
swig_link_libraries(titanlibcustom ${PROJECT_NAME} ${PYTHON_LIBRARIES})
if(APPLE)
    set_property(TARGET _titanlibcustom PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
    set_target_properties(_titanlibcustom PROPERTIES
       SUFFIX ".so"
       INSTALL_RPATH"@loader_path;@loader_path/../../${PROJECT_NAME}/.libs"
       )
endif()
#set_property(TARGET _pytitanlibcustom PROPERTY OUTPUT_NAME titanlibcustom)

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/titanlibcustom.py
        ${CMAKE_CURRENT_BINARY_DIR}/_titanlibcustom.so)

# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

#set(INIT_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/version.in.py)
#set(INIT_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/titanlibcustom/__init__.py)
#configure_file(${INIT_PY_IN} ${INIT_PY_OUT})

# Install target to call setup.py
add_custom_target(install-python
        DEPENDS _titanlibcustom
        COMMAND python3 ${SETUP_PY_OUT} install)

add_custom_target(install-python-user
        DEPENDS _titanlibcustom
        COMMAND python3 ${SETUP_PY_OUT} install --user)

add_custom_target(build-python
        DEPENDS _titanlibcustom
        COMMAND python3 ${SETUP_PY_OUT} build)

add_custom_target(develop-python
        DEPENDS _titanlibcustom
        COMMAND python3 ${SETUP_PY_OUT} develop)

add_custom_target(develop-python-user
        DEPENDS _titanlibcustom
        COMMAND python3 ${SETUP_PY_OUT} develop --user)

