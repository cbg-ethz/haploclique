# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# Adapted by Tobias Marschall when bundling BamTools into HaploClique
#
# src/api/
# ==========================

# list include paths
include_directories( ${BamTools_SOURCE_DIR}/src )

# add compiler definitions 
add_definitions( -DBAMTOOLS_API_LIBRARY ) # (for proper exporting of library symbols)
add_definitions( -fPIC ) # (attempt to force PIC compiling on CentOS, not being set on shared libs by CMake)

# fetch all internal source files
add_subdirectory ( internal )

# make list of all API source files
set( BamToolsAPISources
        BamAlignment.cpp
        BamMultiReader.cpp
        BamReader.cpp
        BamWriter.cpp
        SamHeader.cpp
        SamProgram.cpp
        SamProgramChain.cpp
        SamReadGroup.cpp
        SamReadGroupDictionary.cpp
        SamSequence.cpp
        SamSequenceDictionary.cpp
        ${InternalSources}
)

# create main BamTools API shared library
if(${STATIC_LINKING} MATCHES "1")
	add_library( BamTools STATIC ${BamToolsAPISources} )
else()
	add_library( BamTools SHARED ${BamToolsAPISources} )
endif()
set_target_properties( BamTools PROPERTIES SOVERSION "2.1.0" OUTPUT_NAME "bamtools-haploclique" )

# link libraries automatically with zlib (and Winsock2, if applicable)
if( _WIN32 )
    set( APILibs z ws2_32 )
else( _WIN32 )
    set( APILibs z )
endif( _WIN32 )

target_link_libraries( BamTools ${APILibs} )

# set library install destinations
install( TARGETS BamTools LIBRARY DESTINATION "lib" ARCHIVE DESTINATION "lib")
