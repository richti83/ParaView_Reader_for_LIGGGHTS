#ifndef _liggghts_reader_h
#define _liggghts_reader_h

#include "vtkPolyDataAlgorithm.h"
#include <map> //needed for protected ivars
#include <vector> //needed for protected ivars
#include <string> //needed for protected ivars

#include "vtkInformationIntegerKey.h"

class VTK_EXPORT liggghts_binary_reader : public vtkPolyDataAlgorithm
{
public:
  static liggghts_binary_reader *New();
  vtkTypeMacro(liggghts_binary_reader,vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int CanReadFile(const char* fname);

protected:
  liggghts_binary_reader();
  ~liggghts_binary_reader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int Canreadfile(const char *fname);
  char *FileName;
  char *FieldDelimiterCharacters;
  
  size_t NumberOfPoints;

  void OpenFile();
  ifstream *File;

private:
  liggghts_binary_reader(const liggghts_binary_reader&);
  void operator = (const liggghts_binary_reader&);
};
#endif
