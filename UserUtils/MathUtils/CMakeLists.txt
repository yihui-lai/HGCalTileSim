file(GLOB MathUtils_src src/*.cc)
add_library(MathUtils SHARED ${MathUtils_src})
target_link_libraries( MathUtils
  Common
  ${GSL_LIBRARIES}
  ${ROOT_LIBRARIES}
)

## Function for compiling unit tests
function(make_mathutils_test testfile)
  get_filename_component( testname ${testfile} NAME_WE )
  set( testbin "usrutil_${testname}" )
  add_executable( ${testbin} ${testfile} )
  set_target_properties( ${testbin} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/testbin/MathUtils )
  target_link_libraries( ${testbin} MathUtils )
endfunction()

## Listing all unit test functions to perform
file(GLOB math_test_files "test/*.cc")
foreach(math_test_file ${math_test_files})
  make_mathutils_test( ${math_test_file} )
endforeach()
