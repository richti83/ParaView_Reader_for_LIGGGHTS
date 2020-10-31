set(classes
  liggghts_bondreader)

vtk_module_add_module(liggghts_bondreader_mod
  CLASSES ${classes}
  SOURCES liggghts_bondreader.cxx)

paraview_add_server_manager_xmls(
  XMLS  liggghts_bondreader.xml)