#ifndef _liggghts_reader_h
#define _liggghts_reader_h

#include "vtkPolyDataAlgorithm.h"





#include "vtkInformationIntegerKey.h"

class VTK_EXPORT liggghts_reader : public vtkPolyDataAlgorithm
{
public:
  static liggghts_reader *New();
  vtkTypeMacro(liggghts_reader,vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int CanReadFile(const char* fname);

protected:
  liggghts_reader();
  ~liggghts_reader();
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int Canreadfile(const char *fname);
  char *FileName;
  char *FieldDelimiterCharacters;

  size_t NumberOfPoints;

  void OpenFile();
  ifstream *File;

private:
  liggghts_reader(const liggghts_reader&);
  void operator = (const liggghts_reader&);
};
#endif
