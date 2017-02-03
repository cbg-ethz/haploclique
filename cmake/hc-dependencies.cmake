# External libraries

# Boost
if(NOT Boost_INCLUDE_DIRS OR NOT Boost_LIBRARIES)
    get_property(HC_STATIC GLOBAL PROPERTY HC_STATIC_GLOBAL)
    if (HC_STATIC)
        set(Boost_USE_STATIC_LIBS ON)
    endif()
    find_package(Boost 1.36 COMPONENTS program_options iostreams unit_test_framework REQUIRED)
endif()

# Threads
if (NOT Threads)
    find_package(Threads REQUIRED)
endif()

# ZLIB
if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()