#include "liggghts_binary_reader.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"


#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTable.h"
#include "vtkVariant.h"
#include "vtkObject.h"

#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "vtkObjectFactory.h"
#include "vtkByteSwap.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"
#include "vtkCharArray.h"
#include "vtkFloatArray.h"

#include <algorithm>
#include <vector>
#include <string>
#include <vtksys/ios/sstream>

//#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#ifdef _WIN32
typedef signed __int64       int64_t;
typedef unsigned __int64     uint64_t;
#endif

//vtkCxxRevisionMacro(liggghts_binary_reader, "$Revision: 1.0 $");
vtkStandardNewMacro(liggghts_binary_reader);

liggghts_binary_reader::liggghts_binary_reader()
{
	this->FileName = 0;
	this->File=0;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
}

liggghts_binary_reader::~liggghts_binary_reader()
{
	if (this->File)
	{
		this->File->close();
		delete this->File;
		this->File = NULL;
	}

	this->SetFileName(0);
	this->FileName = NULL;
}

void liggghts_binary_reader::OpenFile()
{
	if (!this->FileName)
	{
		vtkErrorMacro(<<"FileName must be specified.");
		return;
	}
	// If the file was open close it.
	if (this->File)
	{
		this->File->close();
		delete this->File;
		this->File = NULL;
	}

	// Open the new file.

#ifdef _WIN32
	this->File = new ifstream(this->FileName, ios::in | ios::binary);
#else
	this->File = new ifstream(this->FileName, ios::binary);
#endif
	if (! this->File || this->File->fail())
	{
		vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
		return;
	}
}

int liggghts_binary_reader::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	int i,j,m,n;
	int64_t ntimestep, natoms;
	int size_one,nchunk,triclinic;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	int maxbuf = 0;
	double *buf = NULL;

	vtkPoints *points = vtkPoints::New();
	points->SetDataTypeToFloat();
	points->Reset();   

	vtkFloatArray *PointID = vtkFloatArray::New();
	PointID->SetNumberOfComponents(1);
	PointID->Reset();
	PointID->SetName("Point ID");
	PointID->SetComponentName(0,"Point ID");

	vtkFloatArray *v = vtkFloatArray::New();
	v->SetNumberOfComponents(3);
	v->SetComponentName(0,"X");
	v->SetComponentName(1,"Y");
	v->SetComponentName(2,"Z");
	v->SetName("v");

	vtkFloatArray *radius = vtkFloatArray::New();
	radius->SetNumberOfComponents(1);
	radius->SetName("radius");

	vtkFloatArray *scalar = vtkFloatArray::New();
	scalar->SetNumberOfComponents(1);
	scalar->SetName("Scalar");

	double id=0;
	double val[3];
	val[0]=val[1]=val[2]=0;
	double vel[3];
	vel[0]=vel[1]=vel[2]=0;
	double rad=0;
	double scal=0;

	FILE *fp = fopen(this->FileName,"rb");
	if (!fp) {
		vtkErrorMacro(<<"ERROR: Could not open %s\n"<<this->FileName);
		exit(1);
	}   

	while (1) {

		fread(&ntimestep,sizeof(ntimestep),1,fp);

		// detect end-of-file

		if (feof(fp)) {
			fclose(fp);
			break;
		}

		int boundary[6];

		fread(&natoms,sizeof(natoms),1,fp);
		fread(&triclinic,sizeof(int),1,fp);
		fread(&boundary,6*sizeof(int),1,fp);
		fread(&xlo,sizeof(double),1,fp);
		fread(&xhi,sizeof(double),1,fp);
		fread(&ylo,sizeof(double),1,fp);
		fread(&yhi,sizeof(double),1,fp);
		fread(&zlo,sizeof(double),1,fp);
		fread(&zhi,sizeof(double),1,fp);

		fread(&size_one,sizeof(int),1,fp);
		fread(&nchunk,sizeof(int),1,fp);

		//vtkErrorMacro(<<"TIMESTEP" << ntimestep << "#");
		//vtkErrorMacro(<<"NUMBER OF ATOMS" << natoms << "#");
		//vtkErrorMacro(<<"size_one" << size_one << "\n");
		//vtkErrorMacro(<<"nchunk" << nchunk << "\n");

		for (i = 0; i < nchunk; i++) {
			fread(&n,sizeof(int),1,fp); //anzahl gesamt

			// vtkErrorMacro(<<"nchunk" << n << "\n");
			// extend buffer to fit chunk size

			if (n > maxbuf) {
				if (buf) delete [] buf;
				buf = new double[n];
				maxbuf = n;
			}

			// read chunk and write as size_one values per line

			fread(buf,sizeof(double),n,fp); //lies alles 0..n
			n /= size_one; //elemente pro zeile
			m = 0;
			for (j = 0; j < n; j++) {
				//for (k = 0; k < size_one; k++)
				//vtkErrorMacro(<<"j="<<j<<" k="<<k<<" m="<<m<<":"<< buf[m++] << "#");
				//id x y z vx vy vz radius mass
				id=buf[m++];
				PointID->InsertNextTuple1(id);
				val[0]=buf[m++];
				val[1]=buf[m++];
				val[2]=buf[m++];
				points->InsertNextPoint(val[0], val[1], val[2]);
				vel[0]=buf[m++];
				vel[1]=buf[m++];
				vel[2]=buf[m++];
				v->InsertNextTuple(vel);
				radius->InsertNextTuple1(buf[m++]);
				scalar->InsertNextTuple1(buf[m++]);
			}//for j 
		}//for i

	}//while(1)


	vtkCellArray *vertices = vtkCellArray::New();
	vertices->Reset();

	this->NumberOfPoints = points->GetNumberOfPoints();
	for( vtkIdType j = 0; j < (vtkIdType)this->NumberOfPoints; ++j )
	{
		vertices->InsertNextCell( 1 );
		vertices->InsertCellPoint( j );
	}

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the ouptut
	vtkPolyData *myoutput = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	myoutput->SetPoints(points);
	myoutput->SetVerts(vertices);

	myoutput->GetPointData()->AddArray(PointID);
	myoutput->GetPointData()->AddArray(radius);
	myoutput->GetPointData()->AddArray(v);
	myoutput->GetPointData()->AddArray(scalar);


	/////////


	vtkIntArray *intValue = vtkIntArray::New();
	intValue->SetNumberOfComponents(1);
	intValue->SetName("Dumpstep");
	intValue->InsertNextValue(ntimestep);

	vtkIntArray *intValue1 = vtkIntArray::New();
	intValue1->SetNumberOfComponents(1);
	intValue1->SetName("COUNT");
	intValue1->InsertNextValue(natoms);

	vtkStringArray *strValue = vtkStringArray::New();
	strValue->SetNumberOfComponents(1);
	strValue->SetName("fname");
	std::string tstr = this->FileName;
	strValue->InsertNextValue(tstr.c_str());

	myoutput->GetFieldData()->AddArray(intValue);
	myoutput->GetFieldData()->AddArray(intValue1);
	myoutput->GetFieldData()->AddArray(strValue);
	myoutput->Modified();

	// Speicher freigeben
	if (buf)
	{
		delete [] buf;
		buf = NULL;
	}
	PointID->Delete();
	v->Delete();
	radius->Delete();
	scalar->Delete();

	return 1;
}

int liggghts_binary_reader::RequestInformation(
	vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **vtkNotUsed(inputVector),
	vtkInformationVector *outputVector)
{
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	//outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
	outInfo->Set(CAN_HANDLE_PIECE_REQUEST(),
                 1);
	return 1;
}

void liggghts_binary_reader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "FileName: "
		<< (this->FileName ? this->FileName : "(NULL)") << endl;
}

int liggghts_binary_reader::CanReadFile(const char *fname)
{
	return 1;
}
