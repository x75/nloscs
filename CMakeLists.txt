cmake_minimum_required (VERSION 2.8)
set(SUPERNOVA_CMAKE_MINVERSION 3.1)

set(PROJECT "NLOscs")
project (${PROJECT})

set(SC_PATH "$ENV{SC_PATH}")

if(SC_PATH) 
  message("SC_PATH set to ${SC_PATH}")
else()
  message(FATAL_ERROR "Please define the SC_PATH environment variable to point to SuperCollider source folder.")
endif()


include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/common)


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND NOT CMAKE_VERSION VERSION_LESS SUPERNOVA_CMAKE_MINVERSION)
	option(SUPERNOVA "Build plugins for supernova" ON)
else()
	option(SUPERNOVA "Build plugins for supernova" OFF)
endif()

option(NATIVE "Build with native optimization" ON)
option(CPP11 "Build with c++11." ON)
option(OPTIMIZED "Build with -O3" ON)
option(NOVA_SIMD "Build VBAP with nova-simd support." ON)


if (NOVA_SIMD)
  add_definitions(-DNOVA_SIMD)
  include_directories(external_libraries/nova-simd)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_definitions(-fvisibility=hidden)

    include (CheckCCompilerFlag)
    include (CheckCXXCompilerFlag)

    CHECK_C_COMPILER_FLAG(-msse HAS_SSE)
    CHECK_CXX_COMPILER_FLAG(-msse HAS_CXX_SSE)

    if (HAS_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    endif()
    if (HAS_CXX_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
    endif()

    CHECK_C_COMPILER_FLAG(-msse2 HAS_SSE2)
    CHECK_CXX_COMPILER_FLAG(-msse2 HAS_CXX_SSE2)

    if (HAS_SSE2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif()
    if (HAS_CXX_SSE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    endif()

    CHECK_C_COMPILER_FLAG(-mfpmath=sse HAS_FPMATH_SSE)
    CHECK_CXX_COMPILER_FLAG(-mfpmath=sse HAS_CXX_FPMATH_SSE)

    if (HAS_SSE42)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2")
    endif()
    if (HAS_CXX_SSE42)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    endif()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")

    CHECK_C_COMPILER_FLAG(-mfpmath=sse HAS_FPMATH_SSE)
    CHECK_CXX_COMPILER_FLAG(-mfpmath=sse HAS_CXX_FPMATH_SSE)
    if (HAS_FPMATH_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse")
    endif()
    if (HAS_CXX_FPMATH_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
    endif()

    if(NATIVE)
        add_definitions(-march=native)
    endif()

    if(CPP11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        if(CMAKE_COMPILER_IS_CLANG)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif()
      endif()
      
    if(OPTIMIZED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        if(CMAKE_COMPILER_IS_CLANG)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        endif()
    endif()
endif()
if(MINGW)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mstackrealign")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mstackrealign")
endif()

# linear algebra libs

# find_package(Armadillo REQUIRED)
# include_directories(${ARMADILLO_INCLUDE_DIRS})

# flens
# find_path(FLENS_INCLUDE_DIR flens/flens.h ${TRIAL_INCLUDE_PATHS}  NO_DEFAULT_PATH)
# find_path(FLENS_INCLUDE_DIR flens/flens.h)

set(plugin_sources
  NLOscillators.cpp
  )

set(plugins "")
set(supernova_plugins "")

foreach(plugin ${plugin_sources})
  string(REPLACE .cpp "" plugin_name ${plugin} )
  add_library(${plugin_name} MODULE source/${plugin})
  list(APPEND plugins ${plugin_name})
  # target_link_libraries(${plugin_name} ${ARMADILLO_LIBRARIES})
  
  if(SUPERNOVA)
    add_library(${plugin_name}_supernova MODULE source/${plugin})
    # target_link_libraries(${plugin_name}_supernova ${ARMADILLO_LIBRARIES})
    
    list(APPEND supernova_plugins ${plugin_name}_supernova)
    set_property(TARGET ${plugin_name}_supernova
      PROPERTY COMPILE_DEFINITIONS SUPERNOVA)
    if(APPLE)
      install(TARGETS ${plugin_name}_supernova
        DESTINATION ${INSTALL_DESTINATION}/${PLUGIN_DIR})
    endif()
  endif()
endforeach(plugin)



# add_library(EvoPlast MODULE source/EvoPlast.cpp)
# add_library(Reservoir MODULE source/Reservoir.cpp)

# add_executable(foo foo.cc)

