add_subdirectory(model)
add_subdirectory(utils)
pybind11_add_module(python_sbm python_API.cpp)
target_link_libraries(python_sbm PUBLIC
        project_options
        project_warnings
        Boost::system
        Boost::thread
        sbm::model
        sbm::utils
        OpenMP::OpenMP_CXX)

target_include_directories(python_sbm PRIVATE ${PROJECT_SOURCE_DIR}/includes)
target_compile_definitions(python_sbm
        PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

add_executable(whole_blood_sbm main.cpp)
target_link_libraries(whole_blood_sbm PUBLIC
        project_options
        project_warnings
        Boost::system
        Boost::thread
        sbm::model
        sbm::utils
        OpenMP::OpenMP_CXX)

target_include_directories(whole_blood_sbm PRIVATE ${PROJECT_SOURCE_DIR}/includes)


