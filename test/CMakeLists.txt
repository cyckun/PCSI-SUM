add_executable(pcsi_sum
        pcsi_sum_test.cpp
        )

target_link_libraries(pcsi_sum PUBLIC
        pcsi
        )

set_target_properties(pcsi_sum
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )
