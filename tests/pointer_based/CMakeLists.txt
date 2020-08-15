set(test_NAMES
        BlockTreeBasicPropertiesFixture
        BlockTreeFixture
        BlockTreeWithoutCleanningFixture
        HeuristicBlockTreeFixture)


foreach( exe ${test_NAMES} )
    add_executable( check_${exe} ${exe}.cpp)
endforeach()


foreach( exe ${test_NAMES} )
    target_link_libraries( check_${exe} gtest gtest_main)
    target_link_libraries( check_${exe} blocktree )

    target_link_libraries( check_${exe} sdsl )
    target_link_libraries( check_${exe} divsufsort )
    target_link_libraries( check_${exe} divsufsort64 )

    target_link_libraries( check_${exe} gctbio )
    target_link_libraries( check_${exe} textIndexes )
    target_link_libraries( check_${exe} cds )

    target_link_libraries( check_${exe} libRepairGN )
endforeach()