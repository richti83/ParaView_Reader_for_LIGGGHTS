#include "liggghts_bondreader.h"
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

#include <vtkCellData.h>

#include <algorithm>
#include <vector>
#include <string>
#include <vtksys/ios/sstream>

//#include "vtkImageAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"

/*
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
*/



//vtkCxxRevisionMacro(liggghts_bondreader, "$Revision: 2.0 $");
vtkStandardNewMacro(liggghts_bondreader);
//vtkInformationKeyMacro(liggghts_bondreader, TS_KEY, Integer);

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

liggghts_bondreader::liggghts_bondreader()
{
	this->FileName = 0;
	this->File=0;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
}

liggghts_bondreader::~liggghts_bondreader()
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

void liggghts_bondreader::OpenFile()
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

// Entfernt Sonderzeichen am Ende
void trim(char *str)
{
	size_t i = strlen(str)-1;
	while( (i>=0) && ((str[i] == '\r') || (str[i] == '\n') || (str[i] == ' ')))
	{
		str[i] = '\0';
		i--;
	}
}

int liggghts_bondreader::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{

	/*struct timeval begin, end;
	long seconds, useconds;

	gettimeofday(&begin,(struct timezone *)0);
	*/
	float simtime=-1.0;
	this->OpenFile();
	char line[512]; //increase to get more chars !
	std::string label;

	this->File->getline(line,sizeof(line)); //1st line
	this->File->getline(line,sizeof(line)); //2nd line = Timestep
	int TS=atoi(line);

	this->File->getline(line,sizeof(line)); //No of items
	this->File->getline(line,sizeof(line)); //4th line = #Atoms
	int COUNT=atoi(line);

	this->File->getline(line,sizeof(line));
	trim(line);
	if (strncmp(line,"ITEM: BOX BOUNDS",15) == 0) {
		//override ..
		this->File->getline(line,sizeof(line)); //xlow xhi
		this->File->getline(line,sizeof(line)); //ylow yhi
		this->File->getline(line,sizeof(line)); //zlow zhi
		this->File->getline(line,sizeof(line)); //Item Entries
	}	
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetDataTypeToFloat();
	points->Reset();

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

	vtkDoubleArray *pid = vtkDoubleArray::New();
	pid->SetNumberOfComponents(1);
	pid->Reset();
	pid->SetNumberOfTuples(2*COUNT);
	pid->SetComponentName(0,"ID");
	pid->SetName("Particle");

	vtkDoubleArray *mol = vtkDoubleArray::New();
	mol->SetNumberOfComponents(1);
	mol->Reset();
	mol->SetNumberOfTuples(2*COUNT);
	mol->SetName("Molecule");

	vtkDoubleArray *btype = vtkDoubleArray::New();
	btype->SetNumberOfComponents(1);
	btype->Reset();
	btype->SetNumberOfTuples(2*COUNT);
	btype->SetName("Type");

	vtkDoubleArray *force = vtkDoubleArray::New();
	force->SetNumberOfComponents(3);
	force->SetName("F");

	if (COUNT>0) {

	float x1[3],x2[3];
	int btype_,mol_;
	//double fx_,fy_,fz_;
	double f[3];
	int lc=0;
	int pc=0;
	std::string item;
	//printf("expecting %d elements\n",COUNT); 

	// Allocate memory
	int *id1 = new int[COUNT];
	if (id1 == NULL) 
	{
		printf("Error allocating memory for id1!\n"); 
		return 0;
	}
	int *id2 = new int[COUNT];
	if (id2 == NULL) 
	{
		printf("Error allocating memory for id2!\n"); 
		return 0;
	}



	vtkSmartPointer<vtkLine> *vtkline = new vtkSmartPointer<vtkLine>[COUNT];
	if (vtkline == NULL) 
	{
		printf("Error allocating memory for vtkline!\n"); 
		return 0;
	}

	//Create a cell array to store the lines in and add the lines to it


	while ( this->File->getline(line,sizeof(line)) ) { //every line
		int ic=0;
		std::string line_content=line;
		std::stringstream ls(line_content);
		while(std::getline(ls, item, ' ')) { //every item
			switch (ic) {
			case 0:	x1[0]=atof(item.c_str());break;
			case 1:	x1[1]=atof(item.c_str());break;
			case 2:	x1[2]=atof(item.c_str());break;
			case 3:	x2[0]=atof(item.c_str());break;
			case 4:	x2[1]=atof(item.c_str());break;
			case 5:	x2[2]=atof(item.c_str());break;
			case 6: id1[lc]=atoi(item.c_str());break;
			case 7: id2[lc]=atoi(item.c_str());break;
			case 8: btype_=atoi(item.c_str());break;
			case 9: mol_=atoi(item.c_str());break;
			case 10: f[0]=atof(item.c_str());break;
			case 11: f[1]=atof(item.c_str());break;
			case 12: f[2]=atof(item.c_str());break; 
			}
			ic++;
		} //item
		points->InsertNextPoint(x1[0], x1[1], x1[2]);
		points->InsertNextPoint(x2[0], x2[1], x2[2]);

		pid->InsertTuple1(pc, id1[lc]);	   //ID of point1
		pid->InsertTuple1(pc+1, id2[lc]);  //ID of point2

		btype->InsertTuple1(pc, btype_);    //btype of point1
		btype->InsertTuple1(pc+1, btype_);  //btype of point2

		mol->InsertTuple1(pc, mol_);    //mol of point1
		mol->InsertTuple1(pc+1, mol_);  //mol of point2

		vtkline[lc] = vtkSmartPointer<vtkLine>::New();
		vtkline[lc]->GetPointIds()->SetId(0,pc);
		vtkline[lc]->GetPointIds()->SetId(1,pc+1);
		lines->InsertNextCell(vtkline[lc]);

		force->InsertNextTupleValue(f);
		pc+=2;
		lc++;
	}//line


	// free memory
	if (id1 != NULL) 
	{
		delete [] id1;
		id1 = NULL;
	}
	if (id2 != NULL) 
	{
		delete [] id2;
		id2 = NULL;
	}
	if (vtkline != NULL) 
	{
		delete [] vtkline;
		vtkline = NULL;
	}

	} //count>0
	else {
		points->InsertNextPoint(0, 0, 0);
	     }
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the ouptut
	vtkPolyData *myoutput = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));


	myoutput->SetPoints(points);
	if (COUNT>0) {
		myoutput->SetLines(lines);
		//myoutput->SetVerts(vertices);
		myoutput->GetPointData()->AddArray(pid);
		myoutput->GetPointData()->AddArray(btype);
		myoutput->GetPointData()->AddArray(mol);
		myoutput->GetCellData()->AddArray(force);
		myoutput->Modified();
	}


	return 1;
}

int liggghts_bondreader::RequestInformation(
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

void liggghts_bondreader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "FileName: "
		<< (this->FileName ? this->FileName : "(NULL)") << endl;
}

int liggghts_bondreader::CanReadFile(const char *fname)
{
	return 1;
}
