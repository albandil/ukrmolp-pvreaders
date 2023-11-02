cmake_minimum_required(VERSION 3.5)

project(UKRmolpPlugins C CXX)

include(GNUInstallDirs)

find_package(ParaView REQUIRED)

paraview_plugin_scan(
    PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/DipelmReader/paraview.plugin"
    PROVIDES_PLUGINS ukrmolp_readers
    ENABLE_BY_DEFAULT ON
)

paraview_plugin_build(PLUGINS ${ukrmolp_readers})
