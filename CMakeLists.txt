# -*- mode: CMAKE; -*-

cmake_minimum_required(VERSION 2.8)
project(qdigest CXX)
set(CMAKE_CXX_STANDARD 11)

# build version number generator - NICE!
set(QDIGEST_SOURCE
    qdigest.cpp
)

add_library(qdigest SHARED ${QDIGEST_SOURCE})
target_include_directories(qdigest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(qdigeststatic STATIC ${QDIGEST_SOURCE})
target_include_directories(qdigeststatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(qdigeststatic PROPERTIES OUTPUT_NAME qdigest)
