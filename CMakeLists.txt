cmake_minimum_required(VERSION 3.0)
project(Installments_Scheduling)

include_directories(
    ${PROJECT_SOURCE_DIR}/include/
)

add_library(${PROJECT_NAME} SHARED
    models/MISRR.cpp
    models/PMIS.cpp
    models/RRMIS.cpp
    models/Server.cpp
)

add_executable(main main.cpp)

target_link_libraries(main ${PROJECT_NAME})
