cmake_minimum_required(VERSION 2.6)
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED)

project(MOPSI)

add_executable(MOPSI
            main.cpp
            mopsi.cpp mopsi.h)

target_link_libraries(MOPSI /home/quentin/anaconda3/lib/libpython3.5m.a)

ImagineUseModules(MOPSI Graphics)

#cmake_minimum_required(VERSION 3.12)
#file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
#list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
#list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
#find_package(Imagine REQUIRED)

#project(MOPSI)

## get_python_config_var(var_name target_var) executes the Python interpreter
## to retrieve var_name and stores its value to target_var.
#function(get_python_config_var var_name target_var)
#  execute_process(
#    COMMAND python -c
#      "import sysconfig; print(sysconfig.get_config_var(\"${var_name}\"))"
#    OUTPUT_VARIABLE value OUTPUT_STRIP_TRAILING_WHITESPACE)
#  set(${target_var} ${value} PARENT_SCOPE)
#endfunction(get_python_config_var)

#get_python_config_var(INCLDIRSTOMAKE Python_INCLUDE_DIRS)
#separate_arguments(Python_INCLUDE_DIRS)
#get_python_config_var(LIBDIR Python_LIBRARY_DIRS)
#separate_arguments(Python_LIBRARY_DIRS)
#get_python_config_var(BLDLIBRARY Python_LIBRARIES)
#separate_arguments(Python_LIBRARIES)

## link to python by default
#if (Python_LIBRARIES)
#else()
#  set(Python_LIBRARIES "python")
#endif()

## link_directories only works if before add_executable
#link_directories(MOPSI ${Python_LIBRARY_DIRS})
#add_executable(MOPSI
#               main.cpp
#               mopsi.cpp mopsi.h)
#include_directories(MOPSI ${Python_INCLUDE_DIRS})
#target_link_libraries(MOPSI ${Python_LIBRARIES})

#ImagineUseModules(MOPSI Graphics)
