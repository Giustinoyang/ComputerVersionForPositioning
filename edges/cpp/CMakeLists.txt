cmake_minimum_required(VERSION 2.8)

link_directories("/home/xcy/anaconda3/lib")
link_directories("/home/xcy/edges/cpp/build")

add_compile_options( -D_GLIBCXX_USE_CXX11_ABI=0)

IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
    build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
    Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
    options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
    RelWithDebInfo MinSizeRel.")
ENDIF()

project(edge_boxes)

find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(Eigen3)
find_package(Boost COMPONENTS system filesystem python)
find_package(PythonLibs)

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-sign-compare -pedantic") 
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DEIGEN_DONT_PARALLELIZE")

include_directories(${OpenCV_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/include")

add_definitions(-DUSEOMP)
set(CVMATIO_PATH "/home/xcy/cvmatio")
set(NUMPY_OPENCV_CONVERTER_PATH "/home/xcy/numpy-opencv-converter-master")

add_executable(edge_boxes_mat src/edge_boxes_mat.cpp src/edge_nms.cpp src/util.cpp)
target_link_libraries(edge_boxes_mat ${OpenCV_LIBS} ${CVMATIO_PATH}/lib/libcvmatio.so)
target_include_directories(edge_boxes_mat PRIVATE "${CVMATIO_PATH}/include")

add_executable(edge_boxes_demo src/edge_boxes_demo.cpp src/edge_boxes.cpp src/edge_detect.cpp src/edge_nms.cpp src/util.cpp)
target_link_libraries(edge_boxes_demo ${OpenCV_LIBS} ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/contour/libcontour.a ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/imgproc/libimgproc.a ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/learning/liblearning.a ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/util/libutil.a)
target_include_directories(edge_boxes_demo PRIVATE "${PROJECT_SOURCE_DIR}/external/gop_1.3/lib/" ${EIGEN3_INCLUDE_DIR})

if(${Boost_FOUND})
    # Eval 
    add_executable(edge_boxes_eval
        src/edge_boxes_eval.cpp
        src/edge_boxes.cpp
        src/edge_detect.cpp
        src/edge_nms.cpp
        src/util.cpp)
    target_link_libraries(edge_boxes_eval
        ${OpenCV_LIBS}
        ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/contour/libcontour.a
        ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/imgproc/libimgproc.a
        ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/learning/liblearning.a
        ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/util/libutil.a
        ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}
        )
    target_include_directories(edge_boxes_eval PRIVATE
        "${PROJECT_SOURCE_DIR}/external/gop_1.3/lib/"
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        )

    # Python wrapper
    if(${PYTHONLIBS_FOUND})
        add_library(edge_boxes_python SHARED
            src/edge_boxes_python.cpp
            src/edge_boxes.cpp
            src/edge_detect.cpp
            src/edge_nms.cpp
            src/util.cpp)
        set_target_properties(edge_boxes_python
            PROPERTIES
            OUTPUT_NAME edge_boxes_python
            COMPILE_FLAGS "${FASTIDIOUS_FLAGS}"
            LINK_FLAGS -dynamic
            PREFIX ""
            )
        target_link_libraries(edge_boxes_python
            ${OpenCV_LIBS}
            ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/contour/libcontour.a
            ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/imgproc/libimgproc.a
            ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/learning/liblearning.a
            ${PROJECT_SOURCE_DIR}/external/gop_1.3/build/lib/util/libutil.a
            ${Boost_PYTHON_LIBRARY}
            ${PYTHON_LIBRARIES}
            "${NUMPY_OPENCV_CONVERTER_PATH}/build/np_opencv_converter.so"
            )
        target_include_directories(edge_boxes_python PRIVATE
            ${PYTHON_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS}
            ${EIGEN3_INCLUDE_DIR}
            "${PROJECT_SOURCE_DIR}/external/gop_1.3/lib/"
            ${NUMPY_OPENCV_CONVERTER_PATH}
            )
    else(${PYTHONLIBS_FOUND})
        MESSAGE(STATUS "Could not find pythonlibs, not compiling edge_boxes_python")
    endif(${PYTHONLIBS_FOUND})

else(${Boost_FOUND})
    MESSAGE(STATUS "Could not find boost, not compiling edge_boxes_eval")
endif(${Boost_FOUND})
