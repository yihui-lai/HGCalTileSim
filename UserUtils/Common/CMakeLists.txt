file(GLOB Common_src src/*.cc)
add_library(Common SHARED ${Common_src})
target_link_libraries( Common
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES}
  stdc++fs
)

## Function for compiling unit tests
function(make_common_test testfile)
  get_filename_component( testname ${testfile} NAME_WE )
  set( testbin "usrutil_${testname}" )
  add_executable( ${testbin} ${testfile} )
  set_target_properties( ${testbin} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/testbin/Common )
  target_link_libraries( ${testbin} Common )
endfunction()

## Listing all unit test functions to perform
file(GLOB common_test_files "test/*.cc")
foreach(common_test_file ${common_test_files})
  make_common_test( ${common_test_file} )
endforeach()