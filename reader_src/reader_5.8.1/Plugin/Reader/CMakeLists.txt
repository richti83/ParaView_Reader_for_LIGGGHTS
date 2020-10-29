set(classes
  liggghts_reader)

vtk_module_add_module(liggghts_reader_mod
  CLASSES ${classes}
  SOURCES liggghts_reader.cxx)

paraview_add_server_manager_xmls(
  XMLS  liggghts_reader.xml)