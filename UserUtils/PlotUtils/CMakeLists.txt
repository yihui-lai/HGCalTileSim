file(GLOB PlotUtils_src src/*.cc)
add_library(PlotUtils SHARED ${PlotUtils_src})
target_link_libraries(PlotUtils
  Common
  ${ROOT_LIBRARIES}
  ${GHOSTSCRIPT_LIBRARIES}
)

## Function for compiling unit tests
function(make_plotutils_test testfile)
  get_filename_component( testname ${testfile} NAME_WE )
  set( testbin "usrutil_${testname}" )
  add_executable( ${testbin} ${testfile} )
  set_target_properties( ${testbin} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/testbin/PlotUtils )
  target_link_libraries( ${testbin} PlotUtils MathUtils )
endfunction()

## Listing all unit test functions to perform
file(GLOB plot_test_files "test/*.cc")
foreach(plot_test_file ${plot_test_files})
  make_plotutils_test( ${plot_test_file} )
endforeach()

