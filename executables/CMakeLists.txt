set(project_EXECUTABLES
        main)


foreach( exe ${project_EXECUTABLES} )
    add_executable( ${exe} ${exe}.cpp)
endforeach()


foreach( exe ${project_EXECUTABLES} )
    target_link_libraries( ${exe} blocktree )

    target_link_libraries( ${exe} sdsl )
    target_link_libraries( ${exe} divsufsort )
    target_link_libraries( ${exe} divsufsort64 )

    target_link_libraries( ${exe} gctbio )
    target_link_libraries( ${exe} textIndexes )
    target_link_libraries( ${exe} cds )

    target_link_libraries( ${exe} libRepairGN )
endforeach()