add_executable(unit 
	main.cpp
	unit1.cpp
)

IF(${STATIC_LINKING} MATCHES "1")
	target_link_libraries(unit
    libhaplo
    libgtest
    ${Boost_LIBRARIES}
    BamTools-static
)
ELSE()
    	target_link_libraries(unit
    libhaplo
    libgtest
    ${Boost_LIBRARIES}
    BamTools
)
ENDIF()

add_test(NAME unit
         COMMAND unit)
