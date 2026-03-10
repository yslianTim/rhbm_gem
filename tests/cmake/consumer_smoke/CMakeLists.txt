cmake_minimum_required(VERSION 3.18)
project(RHBM_GEM_consumer_smoke LANGUAGES CXX)

find_package(RHBM_GEM CONFIG REQUIRED)

add_executable(rhbm_gem_consumer_smoke main.cpp)
target_link_libraries(rhbm_gem_consumer_smoke PRIVATE RHBM_GEM::rhbm_gem)
