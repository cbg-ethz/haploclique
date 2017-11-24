enable_testing()
add_subdirectory(${HC_EXTERNALDIR}/googletest/googletest external/gtest)

get_property(HC_COMPILE_FLAGS GLOBAL PROPERTY HC_COMPIPLE_FLAGS_GLOBAL)
get_property(HC_LINK_FLAGS GLOBAL PROPERTY HC_LINK_FLAGS_GLOBAL)
get_property(HC_STATIC GLOBAL PROPERTY HC_STATIC_GLOBAL)

file(GLOB HC_TEST_CPP "unit/*.cpp")

add_executable(unit EXCLUDE_FROM_ALL ${HC_TEST_CPP})

target_link_libraries(unit
    hc
    gtest_main
    ${Boost_LIBRARIES}
)
if (HC_STATIC)
    target_link_libraries(unit BamTools-static docopt_s)
else()
    target_link_libraries(unit BamTools docopt_s)
endif()

set_target_properties(unit PROPERTIES COMPILE_FLAGS ${HC_COMPILE_FLAGS})
if (HC_LINK_FLAGS)
    set_target_properties(unit PROPERTIES LINK_FLAGS ${HC_LINK_FLAGS})
endif()

add_custom_target(check
    COMMAND unit
    COMMAND __HC_EXE=${CMAKE_BINARY_DIR}/haploclique cram --xunit-file=${CMAKE_BINARY_DIR}/cram.xml ${HC_TESTDIR}/cram/*.t
    DEPENDS haploclique unit
    WORKING_DIRECTORY ${HC_ROOTDIR})
