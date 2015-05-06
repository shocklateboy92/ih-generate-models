project(ih-generate-models)
cmake_minimum_required(VERSION 2.8.4)

set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "--std=c++14")
set (CMAKE_CXX_COMPILER "clang++")

set(SOURCE_FILES
    main.cpp
    )

set(EXTERNAL_SOURCE_FILES
    lib/pugixml/pugixml.cpp
    lib/fastareader/fastareader.cpp
    lib/fastareader/abstractcharsequence.cpp
    )
add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${EXTERNAL_SOURCE_FILES} common.h model.cpp)

include_directories(${CMAKE_SOURCE_DIR}/lib)
#target_include_directories(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/lib)

