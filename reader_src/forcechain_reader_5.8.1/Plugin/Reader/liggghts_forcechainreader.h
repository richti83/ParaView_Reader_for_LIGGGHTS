#ifndef _liggghts_forcechainreader_h
#define _liggghts_forcechainreader_h

#include "vtkPolyDataAlgorithm.h"
#include <map> //needed for protected ivars
#include <vector> //needed for protected ivars
#include <string> //needed for protected ivars

#include "liggghts_forcechainreader_modModule.h"

#include "vtkInformationIntegerKey.h"

class VTK_EXPORT liggghts_forcechainreader : public vtkPolyDataAlgorithm
{
public:
  static liggghts_forcechainreader *New();
  vtkTypeMacro(liggghts_forcechainreader,vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int CanReadFile(const char* fname);

protected:
  liggghts_forcechainreader();
  ~liggghts_forcechainreader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int Canreadfile(const char *fname);
  char *FileName;
  char *FieldDelimiterCharacters;
  
  size_t NumberOfPoints;

  void OpenFile();
  ifstream *File;

private:
  liggghts_forcechainreader(const liggghts_forcechainreader&);
  void operator = (const liggghts_forcechainreader&);
};
#endif
