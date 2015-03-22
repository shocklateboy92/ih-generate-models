project(ih-generate-models)
cmake_minimum_required(VERSION 2.8)

set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "--std=c++14")

set(SOURCE_FILES
    main.cpp
    )

set(EXTERNAL_SOURCE_FILES
    lib/pugixml/pugixml.cpp
    )
add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${EXTERNAL_SOURCE_FILES})

