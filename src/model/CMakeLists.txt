add_library(model SHARED
        PI_SBM_Solver.cpp)
target_link_libraries(model PRIVATE
        Boost::system
        Boost::thread
        project_options
        project_warnings
        OpenMP::OpenMP_CXX)
add_library(sbm::model ALIAS model)

target_include_directories(model PRIVATE ${PROJECT_SOURCE_DIR}/includes)