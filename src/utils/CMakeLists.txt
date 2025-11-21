add_library(utils SHARED
        io_utils.cpp
)
target_link_libraries(utils PRIVATE
        Boost::system
        Boost::thread
        project_options
        project_warnings
        OpenMP::OpenMP_CXX)
add_library(sbm::utils ALIAS utils)
target_include_directories(utils PRIVATE ${PROJECT_SOURCE_DIR}/includes)