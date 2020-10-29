paraview_add_plugin(liggghts_reader
  VERSION "5.8"
  REQUIRED_ON_SERVER
  REQUIRED_ON_CLIENT
  MODULES liggghts_reader_mod
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")