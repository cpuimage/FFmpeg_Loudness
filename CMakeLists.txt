cmake_minimum_required(VERSION 2.8)
project(Loudness)
SET(CMAKE_BUILD_TYPE "Release")

add_executable(${PROJECT_NAME} main.c loudnorm.c)

target_link_libraries(${PROJECT_NAME} -lm)