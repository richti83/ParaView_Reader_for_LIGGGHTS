set(classes
  liggghts_binary_reader)

vtk_module_add_module(liggghts_binary_reader_mod
  CLASSES ${classes}
  SOURCES liggghts_binary_reader.cxx)

paraview_add_server_manager_xmls(
  XMLS  liggghts_binary_reader.xml)