paraview_add_plugin(DipelmReader
    REQUIRED_ON_SERVER
    REQUIRED_ON_CLIENT
    VERSION "1.0"
    MODULES
        UKRmolp::DipelmReaderVTK
    MODULE_FILES
        vtk/vtk.module
    SERVER_MANAGER_XML
        DipelmReader.xml
)
