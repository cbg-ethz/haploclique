# Linker flags
# if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#     set_property(GLOBAL PROPERTY HC_STATIC_GLOBAL ON)
#     set(HC_LINKER_FLAGS "${HC_LINKER_FLAGS} -static")
#     set(HC_LINKER_FLAGS "${HC_LINKER_FLAGS} -static-libstdc++")
# endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "SunOS")
    add_definitions(-DSUN_OS)
    set(HC_LINKER_FLAGS "-lsocket ${HC_LINKER_FLAGS}")
endif()

set_property(GLOBAL PROPERTY HC_LINK_FLAGS_GLOBAL ${HC_LINKER_FLAGS})

# Compile flags
set_property(GLOBAL PROPERTY HC_COMPIPLE_FLAGS_GLOBAL "-std=c++11 -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas")
