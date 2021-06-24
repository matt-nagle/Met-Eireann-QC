#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")

find_package(Boost  REQUIRED)
find_package(gsl  REQUIRED)

## Add libtitanlibcustom
include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(titanlibcustom ${SOURCES})

# Set up public headers for liba
set_target_properties(titanlibcustom PROPERTIES PUBLIC_HEADER "titanlibcustom.h")

# Create main and link to liba
# add_library(${PROJECT_NAME} titanlibcustom.cpp)
#target_link_libraries(${PROJECT_NAME} titanlibcustom)

# Set up public headers for swig-example
# set_target_properties(titanlibcustom PROPERTIES PUBLIC_HEADER "titanlibcustom.h")

# Add install target
install(TARGETS ${PROJECT_NAME} titanlibcustom
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    PUBLIC_HEADER DESTINATION include)
