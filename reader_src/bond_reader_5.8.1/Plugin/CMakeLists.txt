paraview_add_plugin(liggghts_bondreader
  VERSION "5.8"
  REQUIRED_ON_SERVER
  REQUIRED_ON_CLIENT
  MODULES liggghts_bondreader_mod
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")