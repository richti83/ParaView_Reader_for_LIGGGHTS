#include "liggghts_rigids_reader.h"
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
#include "vtkLine.h"

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

#include <vtkCellData.h>

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

//#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

/*
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
*/



//vtkCxxRevisionMacro(liggghts_rigids_reader, "$Revision: 2.0 $");
vtkStandardNewMacro(liggghts_rigids_reader);
//vtkInformationKeyMacro(liggghts_rigids_reader, TS_KEY, Integer);

//tiny little helper
void searchAndReplace(std::string& value, std::string const& search,std::string const& replace)
{
	std::string::size_type  next;

	for(next = value.find(search);        // Try and find the first match
		next != std::string::npos;        // next is npos if nothing was found
		next = value.find(search,next)    // search for the next match starting after
		// the last match that was found.
		)
	{
		// Inside the loop. So we found a match.
		value.replace(next,search.length(),replace);   // Do the replacement.
		next += replace.length();                      // Move to just after the replace
		// This is the point were we start
		// the next search from. 
	}
}

liggghts_rigids_reader::liggghts_rigids_reader()
{
	this->FileName = 0;
	this->File=0;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
}

liggghts_rigids_reader::~liggghts_rigids_reader()
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

void liggghts_rigids_reader::OpenFile()
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
	this->File = new ifstream(this->FileName, ios::in);
#endif
	if (! this->File || this->File->fail())
	{
		vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
		return;
	}
}

int liggghts_rigids_reader::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{

	vtkPoints *points;
	vtkCellArray *vertices;

	points = vtkPoints::New();
	points->SetDataTypeToFloat();
	points->Reset();

	this->OpenFile();
	char line[512]; //increase to get more chars !
	std::string label;

	this->File->getline(line,sizeof(line)); //1st line
	this->File->getline(line,sizeof(line)); //2nd line = Timestep
	int TS=atoi(line);

	this->File->getline(line,sizeof(line)); //No of items
	this->File->getline(line,sizeof(line)); //4th line = #Atoms
	int COUNT=atoi(line);
	//if (COUNT<1) return 1;
	this->File->getline(line,sizeof(line)); //5th line = ITEM: BOX Bounds OR ITEM: ENTRIES
	//trim(line);
	if (strncmp(line,"ITEM: BOX BOUNDS",15) == 0) {
		//override ..
		this->File->getline(line,sizeof(line)); //xlow xhi
		this->File->getline(line,sizeof(line)); //ylow yhi
		this->File->getline(line,sizeof(line)); //zlow zhi
		this->File->getline(line,sizeof(line)); //Item Entries
	}


	double x[3],F[3],Q[4],vel[3],M[9],rot[3],axis[3],angle;
	int id,type;
	double sc;
	int lc=0;

	std::string item;
	//printf("expecting %d elements\n",COUNT); 

	// Setup the point arrays
	vtkIntArray *ids = vtkIntArray::New();
	ids->SetNumberOfComponents(1);
	ids->SetName("id");

	vtkDoubleArray *scalars = vtkDoubleArray::New();
	scalars->SetNumberOfComponents(1);
	scalars->SetName("sc");

	vtkIntArray *types = vtkIntArray::New();
	types->SetNumberOfComponents(1);
	types->SetName("type");

	vtkDoubleArray *forces = vtkDoubleArray::New();
	forces->SetNumberOfComponents(3);
	forces->SetName("F");

	vtkDoubleArray *velocity = vtkDoubleArray::New();
	velocity->SetNumberOfComponents(3);
	velocity->SetName("v");

	vtkDoubleArray *quats = vtkDoubleArray::New();
	quats->SetNumberOfComponents(4);
	quats->SetName("Q");

	vtkDoubleArray *T = vtkDoubleArray::New();
	T->SetNumberOfComponents(9);
	T->SetName("M");

	vtkDoubleArray *R = vtkDoubleArray::New();
	R->SetNumberOfComponents(3);
	R->SetName("ROT");

	vtkDoubleArray *Roll = vtkDoubleArray::New();
	Roll->SetNumberOfComponents(1);
	Roll->SetName("Roll");

	vtkDoubleArray *Pitch = vtkDoubleArray::New();
	Pitch->SetNumberOfComponents(1);
	Pitch->SetName("Pitch");

	vtkDoubleArray *Yaw = vtkDoubleArray::New();
	Yaw->SetNumberOfComponents(1);
	Yaw->SetName("Yaw");

	vtkDoubleArray *A = vtkDoubleArray::New();
	A->SetNumberOfComponents(3);
	A->SetName("AXIS");

	vtkDoubleArray *a = vtkDoubleArray::New();
	a->SetNumberOfComponents(1);
	a->SetName("angle");

        //c_bid c_btype c_xcm[1] c_xcm[2] c_xcm[3] c_quat[1] c_quat[2] c_quat[3] c_quat[4] c_vel[1] c_vel[2] c_vel[3] c_fcm[1] c_fcm[2] c_fcm[3] 
	while ( this->File->getline(line,sizeof(line)) ) { //every line
		int ic=0;
		std::string line_content=line;
		std::stringstream ls(line_content);
		while(std::getline(ls, item, ' ')) { //every item
			switch (ic) {
			case 0:	id=atoi(item.c_str());break;
			case 1:	type=atoi(item.c_str());break;
			case 2:	x[0]=atof(item.c_str());break;
			case 3:	x[1]=atof(item.c_str());break;
			case 4:	x[2]=atof(item.c_str());break;
			case 5:	Q[0]=atof(item.c_str());break;
			case 6:	Q[1]=atof(item.c_str());break;
			case 7:	Q[2]=atof(item.c_str());break;
			case 8:	Q[3]=atof(item.c_str());break;
			case 9:	 vel[0]=atof(item.c_str());break;
			case 10: vel[1]=atof(item.c_str());break;
			case 11: vel[2]=atof(item.c_str());break;
			case 12: F[0]=atof(item.c_str());break;
			case 13: F[1]=atof(item.c_str());break;
			case 14: F[2]=atof(item.c_str());break;
			case 15: sc=atof(item.c_str());break;
			}
			ic++;
		} //item
		points->InsertNextPoint(x[0], x[1], x[2]);
		//printf("lc=%d:%s\n",lc,line);
		forces->InsertNextTupleValue(F);
		quats->InsertNextTupleValue(Q);
		double x,y,z;

		double s=sqrt(1-Q[0]*Q[0]);
		axis[0]=s>0.001?Q[1]/s:1;
		axis[1]=s>0.001?Q[2]/s:0;
		axis[2]=s>0.001?Q[3]/s:0;
		A->InsertNextTupleValue(axis);

		double angle=2*acos(Q[0]);
		//if (angle>=M_PI) angle=M_PI-angle-M_PI/100;
		a->InsertNextValue(angle);
		//cos(a/2) + i ( x * sin(a/2)) + j (y * sin(a/2)) + k ( z * sin(a/2))
		/*Q[0]=cos(angle/2);
		Q[1]=axis[0]*sin(angle/2);
		Q[2]=axis[1]*sin(angle/2);
		Q[3]=axis[2]*sin(angle/2);*/
				

		double wq = Q[0];
		double xq = Q[1];
		double yq = Q[2];
		double zq = Q[3];

		double q0 = Q[0];
		double q1 = Q[1];
		double q2 = Q[2];
		double q3 = Q[3];


		M[0]=xq*xq+wq*wq-yq*yq-zq*zq;
		M[1]=2*(xq*yq-wq*zq);
                M[2]=2*(xq*zq+wq*yq);
                M[3]=2*(wq*zq+xq*yq);
		M[4]=wq*wq-xq*xq+yq*yq-zq*zq;
		M[5]=2*(yq*zq-wq*xq);
		M[6]=2*(xq*zq-wq*yq);
		M[7]=2*(wq*xq+yq*zq);
		M[8]=wq*wq-xq*xq-yq*yq+zq*zq;
		T->InsertNextTupleValue(M);

	        //(Q_1^2+Q_0^2-Q_2^2-Q_3^3)*iHat+2*(Q_0*Q_3+Q_1*Q_2)*jHat+2*(Q_1*Q_3-Q_0*Q_2)*kHat
				
		rot[0]=atan2(2*(q0*q1+q2*q3),1-2*(q1*q1+q2*q2)); //M[0];
		rot[1]=asin(2*(q0*q2-q3*q1));//M[3];
		rot[2]=atan2(2*(q0*q3+q1*q2),1-2*(q2*q2+q3*q3));//M[6];
		R->InsertNextTupleValue(rot);

		double roll_value=atan2(2*yq*wq - 2*xq*zq, 1 - 2*yq*yq - 2*zq*zq);
		double pitch_value=atan2(2*xq*wq - 2*yq*zq, 1 - 2*xq*xq - 2*zq*zq);
		double yaw_value=asin(2*xq*yq + 2*zq*wq);

		Roll->InsertNextValue(roll_value);
		Pitch->InsertNextValue(pitch_value);
		Yaw->InsertNextValue(yaw_value);
		
		ids->InsertNextValue(id);
		scalars->InsertNextValue(sc);
		types->InsertNextValue(type);
		velocity->InsertNextTupleValue(vel);
		lc++;
	}//line

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the ouptut
	vtkPolyData *myoutput = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));


	vertices = vtkCellArray::New();
	vertices->Reset();
        this->NumberOfPoints = points->GetNumberOfPoints();
        //printf("%d points\n",this->NumberOfPoints);
	for( vtkIdType j = 0; j < (vtkIdType)this->NumberOfPoints; ++j )
	    {
		    vertices->InsertNextCell( 1 );
		    vertices->InsertCellPoint( j );
	    }

	myoutput->SetPoints(points);
  	myoutput->SetVerts(vertices);
	myoutput->GetPointData()->AddArray(ids);
	myoutput->GetPointData()->AddArray(scalars);
	myoutput->GetPointData()->AddArray(types);
	myoutput->GetPointData()->AddArray(forces);
	myoutput->GetPointData()->AddArray(quats);
	myoutput->GetPointData()->AddArray(Roll);
	myoutput->GetPointData()->AddArray(Pitch);
	myoutput->GetPointData()->AddArray(Yaw);
	myoutput->GetPointData()->AddArray(velocity);
	myoutput->GetPointData()->AddArray(T);
	myoutput->GetPointData()->AddArray(R);
	myoutput->GetPointData()->AddArray(A);
	myoutput->GetPointData()->AddArray(a);
	
	// free memory
	points->Delete();
	vertices->Delete();
	forces->Delete();forces=NULL;
	quats->Delete();quats=NULL;
	velocity->Delete();velocity=NULL;
	T->Delete();T=NULL;
	R->Delete();R=NULL;
	A->Delete();A=NULL;
	a->Delete();a=NULL;
	ids->Delete();ids=NULL;
	scalars->Delete();scalars=NULL;
	types->Delete();types=NULL;
	Roll->Delete();Roll=NULL;
	Pitch->Delete();Pitch=NULL;
	Yaw->Delete();Yaw=NULL;

	vtkIntArray *intValue;
	intValue=vtkIntArray::New();
	intValue->SetNumberOfComponents(1);
	intValue->SetName("Dumpstep");
	intValue->InsertNextValue(TS);
	myoutput->GetFieldData()->AddArray(intValue);
	intValue->Delete();

	intValue = vtkIntArray::New();
	intValue->SetNumberOfComponents(1);
	intValue->SetName("COUNT");
	intValue->InsertNextValue(COUNT);
	myoutput->GetFieldData()->AddArray(intValue);
	intValue->Delete();


 return 1;
}

int liggghts_rigids_reader::RequestInformation(
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

void liggghts_rigids_reader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "FileName: "
		<< (this->FileName ? this->FileName : "(NULL)") << endl;
}

int liggghts_rigids_reader::CanReadFile(const char *fname)
{
	return 1;
}
