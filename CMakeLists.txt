cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

rosbuild_init()
#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#rosbuild_genmsg()
#rosbuild_gensrv()

include_directories("/usr/include/python2.7")


include_directories(include)
link_directories(${PROJECT_SOURCE_DIR}/lib)

add_executable(driver src/driver.cpp)
target_link_libraries(driver glut sixense_x64 sixense_utils_x64 ${ROS_LIBRARIES} boost_python python2.7)

##add_executable(dlopen_test src/dlopen_test.cpp)
#target_link_libraries(dlopen_test dl)

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})
