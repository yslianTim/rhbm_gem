file(GLOB CONFIGURE_DEPENDS SOURCE *.cpp)

foreach(test_src ${SOURCE})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_project_test(${test_name} ${test_src})
endforeach()
