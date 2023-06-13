cmake_minimum_required(VERSION 3.0)
project(Installments_Scheduling)

include_directories(
    ${PROJECT_SOURCE_DIR}/include/
)

add_library(${PROJECT_NAME} SHARED
    models/MISRRLL.cpp
    models/MISRRL.cpp
    models/myAPMISRR.cpp
    models/APMISRR.cpp
    models/MISRR.cpp
    models/SIS.cpp
    models/PMIS.cpp
    models/RRMIS.cpp
    models/Server.cpp
    models/method.cpp
    exp/exp_1.cpp
)

add_executable(main main.cpp)

target_link_libraries(main ${PROJECT_NAME})
