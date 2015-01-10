#ifndef _liggghts_bondreader_h
#define _liggghts_bondreader_h

#include "vtkPolyDataAlgorithm.h"
#include <map> //needed for protected ivars
#include <vector> //needed for protected ivars
#include <string> //needed for protected ivars

#include "vtkInformationIntegerKey.h"

class VTK_EXPORT liggghts_bondreader : public vtkPolyDataAlgorithm
{
public:
  static liggghts_bondreader *New();
  vtkTypeMacro(liggghts_bondreader,vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int CanReadFile(const char* fname);

protected:
  liggghts_bondreader();
  ~liggghts_bondreader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int Canreadfile(const char *fname);
  char *FileName;
  char *FieldDelimiterCharacters;
  
  size_t NumberOfPoints;

  void OpenFile();
  ifstream *File;

private:
  liggghts_bondreader(const liggghts_bondreader&);
  void operator = (const liggghts_bondreader&);
};
#endif
