set(classes
  liggghts_forcechainreader)

vtk_module_add_module(liggghts_forcechainreader_mod
  CLASSES ${classes}
  SOURCES liggghts_forcechainreader.cxx)

paraview_add_server_manager_xmls(
  XMLS liggghts_forcechainreader.xml)