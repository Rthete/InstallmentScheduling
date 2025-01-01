cmake_minimum_required(VERSION 3.8)
project(Installments_Scheduling)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    ${PROJECT_SOURCE_DIR}/include/
)

include_directories(/usr/include/eigen3)

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
    exp/exp_2.cpp
    exp/exp_3.cpp
    exp/exp_4.cpp
)

add_executable(main main.cpp)

target_link_libraries(main ${PROJECT_NAME})
