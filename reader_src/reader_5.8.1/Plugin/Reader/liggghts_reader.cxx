//VERSION vom 01-02-2014 mit verbesserter Speicherverwaltung OHNE lump enumeration
//#define MEASURE_TIME

#include "liggghts_reader.h"
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
//#include "vtkSmartPointer.h"
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

#include <sstream>
#include "vtkPolyDataAlgorithm.h"

#if defined(WIN32) && defined(MEASURE_TIME)
#include <Windows.h>
#endif

#define MAX(A,B) ((A) > (B) ? (A) : (B))

//vtkCxxRevisionMacro(liggghts_reader, "$Revision: 2.0 $");
vtkStandardNewMacro(liggghts_reader);
//vtkInformationKeyMacro(liggghts_reader, TS_KEY, Integer);

//tiny little helper
void searchAndReplace(std::string& value, std::string const& search, std::string const& replace)
{
	std::string::size_type  next;

	for (next = value.find(search);        // Try and find the first match
		next != std::string::npos;        // next is npos if nothing was found
		next = value.find(search, next)    // search for the next match starting after
										   // the last match that was found.
		)
	{
		// Inside the loop. So we found a match.
		value.replace(next, search.length(), replace);   // Do the replacement.
		next += replace.length();                      // Move to just after the replace
													   // This is the point were we start
													   // the next search from.
	}
}

liggghts_reader::liggghts_reader()
{
	this->FileName = 0;
	this->File = 0;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
}

liggghts_reader::~liggghts_reader()
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

void liggghts_reader::OpenFile()
{
	if (!this->FileName)
	{
		vtkErrorMacro(<< "FileName must be specified.");
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
	if (!this->File || this->File->fail())
	{
		vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
		return;
	}
}

// Entfernt Sonderzeichen am Ende
void trim(char *str)
{
	size_t i = strlen(str) - 1;
	while ((i >= 0) && ((str[i] == '\r') || (str[i] == '\n') || (str[i] == ' ')))
	{
		str[i] = '\0';
		i--;
	}
}

int liggghts_reader::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	vtkPoints *points = NULL;
	vtkCellArray *vertices = NULL;

	vtkFloatArray **vectors = NULL;
	vtkFloatArray **scalars = NULL;
	vtkFloatArray **others = NULL;

#define MSG(x) vtkWarningMacro(x)

#if defined(WIN32) && defined(MEASURE_TIME)
	LARGE_INTEGER start, stop, LI_freq;
	double freq;
	QueryPerformanceFrequency(&LI_freq);
	freq = (double)LI_freq.QuadPart;
	QueryPerformanceCounter(&start);
#endif

	float simtime = -1.0;
	this->OpenFile();
	char line[16384]; //increase to get more chars !
	std::string label;

	this->File->getline(line, sizeof(line)); //1st line
	this->File->getline(line, sizeof(line)); //2nd line = Timestep
	int TS = atoi(line);

	this->File->getline(line, sizeof(line)); //3rd line could be SIMTIME or Numb. of. Atoms

											 // Etnfernen von Sonderzeichen am Ende von Line, ansonsten kann es sein, dass strcmp nicht funktioniert
	trim(line);

	if (strcmp(line, "ITEM: LABEL") == 0) { //backwardscompatibility to old dumpfiles
		this->File->getline(line, sizeof(line)); //4th line = SIMTIME if (a)
		label = line;
		this->File->getline(line, sizeof(line)); //5th line = Headertext Atoms
	}
	this->File->getline(line, sizeof(line)); //3rd || 5th line = #Atoms
	int COUNT = atoi(line);

	for (int i = 0; i < 4; i++) { //override next 4 lines (Box Bounds) ..
		this->File->getline(line, sizeof(line));
	}

	this->File->getline(line, sizeof(line)); //11th line = content
	std::string content = line;
	content = content.substr(12, content.length()); //cut first 12 chars (ITEM: ATOMS)

	std::stringstream ss(content);
	std::string item;
	int pos = 0;
	int nov = 0; //number of vectors
	int nos = 0; //number of scalars
	int noo = 0; //number of other variants (fixes, computes, variables, colors ...)

	int xpos = -1;
	int ypos = -1;
	int zpos = -1;
	int IDpos = -1;
	int MOLpos = -1;
	int TYPEpos = -1;
	int MASSpos = -1;
	int RADpos = -1;
	int Qpos = -1;
	int DENSpos = -1;
	int NUMBONDpos = -1;

	int vpos[3] = { -1,-1,-1 };
	int fpos[3] = { -1,-1,-1 };
	int ipos[3] = { -1,-1,-1 };
	int opos[3] = { -1,-1,-1 };

	int XSpos[3] = { -1,-1,-1 };
	int XUpos[3] = { -1,-1,-1 };
	int MUpos[3] = { -1,-1,-1 };

	int ANGMOMpos[3] = { -1,-1,-1 };
	int QUATpos[4] = { -1,-1,-1,-1 };
	int TQpos[3] = { -1,-1,-1 };
	int OTHERpos[128] = { -1 };

	std::string OTHERname[128];

	while (std::getline(ss, item, ' ')) {
		if (item.compare("id") == 0) {
			//vtkErrorMacro(<<"found ID@POS" << pos);
			IDpos = pos;
			nos++;
			pos++;
		}
		else if (item.compare("mol") == 0) {
			//vtkErrorMacro(<<"found mol@POS" << pos);
			MOLpos = pos;
			nos++;
			pos++;
		}
		else if (item.compare("type") == 0) {
			//vtkErrorMacro(<<"found type@POS" << pos);
			TYPEpos = pos;
			nos++;
			pos++;
		}
		else if (item.compare("mass") == 0) {
			//vtkErrorMacro(<<"found mass@POS" << pos);
			MASSpos = pos;
			nos++;
			pos++;
		}
		else if (item.compare("x") == 0) {
			//vtkErrorMacro(<<"found x@POS" << pos);
			xpos = pos;
			pos++;
		}
		else if (item.compare("y") == 0) {
			//vtkErrorMacro(<<"found y@POS" << pos);
			ypos = pos;
			pos++;
		}
		else if (item.compare("z") == 0) {
			//vtkErrorMacro(<<"found z@POS" << pos);
			zpos = pos;
			pos++;
		}
		else if (item.compare("ix") == 0) {
			//vtkErrorMacro(<<"found ix@POS" << pos);
			nov++;
			ipos[0] = pos;
			pos++;
		}
		else if (item.compare("iy") == 0) {
			//vtkErrorMacro(<<"found iy@POS" << pos);
			ipos[1] = pos;
			pos++;
		}
		else if (item.compare("iz") == 0) {
			//vtkErrorMacro(<<"found iz@POS" << pos);
			ipos[2] = pos;
			pos++;
		}
		else if (item.compare("vx") == 0) {
			//vtkErrorMacro(<<"found vx@POS" << pos);
			nov++;
			vpos[0] = pos;
			pos++;
		}
		else if (item.compare("vy") == 0) {
			//vtkErrorMacro(<<"found vy@POS" << pos);
			vpos[1] = pos;
			pos++;
		}
		else if (item.compare("vz") == 0) {
			//vtkErrorMacro(<<"found vz@POS" << pos);
			vpos[2] = pos;
			pos++;

		}
		else if (item.compare("xs") == 0) {
			//vtkErrorMacro(<<"found xs@POS" << pos);
			nov++;
			XSpos[0] = pos;
			pos++;
		}
		else if (item.compare("ys") == 0) {
			//vtkErrorMacro(<<"found ys@POS" << pos);
			XSpos[1] = pos;
			pos++;
		}
		else if (item.compare("zs") == 0) {
			//vtkErrorMacro(<<"found zs@POS" << pos);
			XSpos[2] = pos;
			pos++;

		}
		else if (item.compare("xu") == 0) {
			//vtkErrorMacro(<<"found xu@POS" << pos);
			nov++;
			XUpos[0] = pos;
			pos++;
		}
		else if (item.compare("yu") == 0) {
			//vtkErrorMacro(<<"found yu@POS" << pos);
			XUpos[1] = pos;
			pos++;
		}
		else if (item.compare("zu") == 0) {
			//vtkErrorMacro(<<"found zu@POS" << pos);
			XUpos[2] = pos;
			pos++;

		}
		else if (item.compare("mux") == 0) {
			//vtkErrorMacro(<<"found mux@POS" << pos);
			nov++;
			MUpos[0] = pos;
			pos++;
		}
		else if (item.compare("muy") == 0) {
			//vtkErrorMacro(<<"found muy@POS" << pos);
			MUpos[1] = pos;
			pos++;
		}
		else if (item.compare("muz") == 0) {
			//vtkErrorMacro(<<"found muz@POS" << pos);
			MUpos[2] = pos;
			pos++;

		}
		else if (item.compare("angmomx") == 0) {
			//vtkErrorMacro(<<"found angmomx@POS" << pos);
			nov++;
			ANGMOMpos[0] = pos;
			pos++;
		}
		else if (item.compare("angmomy") == 0) {
			//vtkErrorMacro(<<"found angmomy@POS" << pos);
			ANGMOMpos[1] = pos;
			pos++;
		}
		else if (item.compare("angmomz") == 0) {
			//vtkErrorMacro(<<"found angmomz@POS" << pos);
			ANGMOMpos[2] = pos;
			pos++;

		}
		else if (item.compare("quatw") == 0) {
			//vtkErrorMacro(<<"found quatw@POS" << pos);
			nov++;
			QUATpos[0] = pos;
			pos++;
		}
		else if (item.compare("quati") == 0) {
			//vtkErrorMacro(<<"found quati@POS" << pos);
			QUATpos[1] = pos;
			pos++;
		}
		else if (item.compare("quatj") == 0) {
			//vtkErrorMacro(<<"found quatj@POS" << pos);
			QUATpos[2] = pos;
			pos++;
		}
		else if (item.compare("quatk") == 0) {
			//vtkErrorMacro(<<"found quatk@POS" << pos);
			QUATpos[3] = pos;
			pos++;

		}
		else if (item.compare("tqx") == 0) {
			//vtkErrorMacro(<<"found tqx@POS" << pos);
			nov++;
			TQpos[0] = pos;
			pos++;
		}
		else if (item.compare("tqy") == 0) {
			//vtkErrorMacro(<<"found tqy@POS" << pos);
			TQpos[1] = pos;
			pos++;
		}
		else if (item.compare("tqz") == 0) {
			//vtkErrorMacro(<<"found tqz@POS" << pos);
			TQpos[2] = pos;
			pos++;

		}
		else if (item.compare("fx") == 0) {
			//vtkErrorMacro(<<"found fx@POS" << pos);
			nov++;
			fpos[0] = pos;
			pos++;
		}
		else if (item.compare("fy") == 0) {
			//vtkErrorMacro(<<"found fy@POS" << pos);
			fpos[1] = pos;
			pos++;
		}
		else if (item.compare("fz") == 0) {
			//vtkErrorMacro(<<"found fz@POS" << pos);
			fpos[2] = pos;
			pos++;
		}
		else if (item.compare("omegax") == 0) {
			//vtkErrorMacro(<<"found ox@POS" << pos);
			nov++;
			opos[0] = pos;
			pos++;
		}
		else if (item.compare("omegay") == 0) {
			//vtkErrorMacro(<<"found oy@POS" << pos);
			opos[1] = pos;
			pos++;
		}
		else if (item.compare("omegaz") == 0) {
			//vtkErrorMacro(<<"found oz@POS" << pos);
			opos[2] = pos;
			pos++;
		}
		else if (item.compare("radius") == 0) {
			//vtkErrorMacro(<<"found radius@POS" << pos);
			nos++;
			RADpos = pos;
			pos++;
		}
		else if (item.compare("q") == 0) {
			//vtkErrorMacro(<<"found q@POS" << pos);
			nos++;
			Qpos = pos;
			pos++;
		}
		else if (item.compare("density") == 0) {
			//vtkErrorMacro(<<"found density@POS" << pos);
			nos++;
			DENSpos = pos;
			pos++;
		}
		else if (item.compare("numbonds") == 0) {
			//vtkErrorMacro(<<"found density@POS" << pos);
			nos++;
			NUMBONDpos = pos;
			pos++;
		}
		else {
			//vtkErrorMacro(<<"unknown item @Pos" << pos << " " <<item.c_str());
			OTHERpos[noo] = pos;
			OTHERname[noo] = item.c_str();
			noo++;
			pos++;
		}
	}

	//check if we have 3D coords ..
	if ((xpos<0) && (ypos<0) && (zpos<0)) {
		vtkErrorMacro("FATAL: missing 3D Positions !");
		return 0;
	}

	/*vtkErrorMacro(<<"got "<<nos<<" Scalars");
	vtkErrorMacro(<<"got "<<nov<<" Vecors");
	vtkErrorMacro(<<"got "<<nob<<" Bond Partners");*/

	points = vtkPoints::New();
	points->SetDataTypeToFloat();
	points->Reset();
	points->SetNumberOfPoints(COUNT);

	if (nov>0)
	{
		vectors = new vtkFloatArray*[nov];
		if (vectors == NULL) vtkErrorMacro("FATAL: Not enough Memory for vectors!");
	}
	if (nos>0)
	{
		scalars = new vtkFloatArray*[nos];
		if (scalars == NULL) vtkErrorMacro("FATAL: Not enough Memory for scalars!");
	}
	if (noo > 0)
	{
		others = new vtkFloatArray*[noo];
		if (others == NULL) vtkErrorMacro("FATAL: Not enough Memory for others!");
	}

	int n_vec = 0;
	int n_scal = 0;
	int n_other = 0;

	int s_IDpos = -1;
	int s_TYPEpos = -1;
	int s_RADpos = -1;
	int s_MOLpos = -1;
	int s_MASSpos = -1;
	int s_Qpos = -1;
	int s_DENSpos = -1;
	int s_NUMBONDpos = -1;


	int v_vpos = -1;
	int v_fpos = -1;
	int v_ipos = -1;
	int v_opos = -1;

	int v_XSpos = -1;
	int v_XUpos = -1;
	int v_MUpos = -1;
	int v_ANGMOMpos = -1;
	int v_QUATpos = -1;
	int v_TQpos = -1;

	int s_other[128] = { -1 };

	if (IDpos >= 0) {
		//vtkErrorMacro(<<"ID");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("LMP ID");
		scalars[n_scal]->SetComponentName(0, "LMP ID");
		s_IDpos = n_scal;
		n_scal++;
	}

	if (MOLpos >= 0) {
		//vtkErrorMacro(<<"MOL");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("Mol");
		scalars[n_scal]->SetComponentName(0, "Mol");
		s_MOLpos = n_scal;
		n_scal++;
	}

	if (MASSpos >= 0) {
		//vtkErrorMacro(<<"MASS");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("Mass");
		scalars[n_scal]->SetComponentName(0, "Mass");
		s_MASSpos = n_scal;
		n_scal++;
	}

	if (Qpos >= 0) {
		//vtkErrorMacro(<<"Q");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("Q");
		scalars[n_scal]->SetComponentName(0, "Q");
		s_Qpos = n_scal;
		n_scal++;
	}

	if (DENSpos >= 0) {
		//vtkErrorMacro(<<"DENS");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("DENSITY");
		scalars[n_scal]->SetComponentName(0, "DENSITY");
		s_DENSpos = n_scal;
		n_scal++;
	}

	if (TYPEpos >= 0) {
		//vtkErrorMacro(<<"Type");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("Type");
		scalars[n_scal]->SetComponentName(0, "Type");
		s_TYPEpos = n_scal;
		n_scal++;
	}

	if ((ipos[0] >= 0) && (ipos[1] >= 0) && (ipos[2] >= 0)) {
		//vtkErrorMacro(<<"I"<<ipos[0]<<ipos[1]<<ipos[2]);
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("I");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_ipos = n_vec;
		n_vec++;
	}

	if ((XSpos[0] >= 0) && (XSpos[1] >= 0) && (XSpos[2] >= 0)) {
		//vtkErrorMacro(<<"xs");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "xs");
		vectors[n_vec]->SetComponentName(1, "ys");
		vectors[n_vec]->SetComponentName(2, "zs");
		vectors[n_vec]->SetName("S");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_XSpos = n_vec;
		n_vec++;
	}

	if ((XUpos[0] >= 0) && (XUpos[1] >= 0) && (XUpos[2] >= 0)) {
		//vtkErrorMacro(<<"xu");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "xu");
		vectors[n_vec]->SetComponentName(1, "yu");
		vectors[n_vec]->SetComponentName(2, "zu");
		vectors[n_vec]->SetName("U");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_XUpos = n_vec;
		n_vec++;
	}

	if ((MUpos[0] >= 0) && (MUpos[1] >= 0) && (MUpos[2] >= 0)) {
		//vtkErrorMacro(<<"mu");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("mu");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_MUpos = n_vec;
		n_vec++;
	}

	if ((ANGMOMpos[0] >= 0) && (ANGMOMpos[1] >= 0) && (ANGMOMpos[2] >= 0)) {
		//vtkErrorMacro(<<"ANGMOM");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("angmom");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_ANGMOMpos = n_vec;
		n_vec++;
	}

	if ((QUATpos[0] >= 0) && (QUATpos[1] >= 0) && (QUATpos[2] >= 0) && (QUATpos[3] >= 0)) {
		//vtkErrorMacro(<<"QUAT");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(4);
		vectors[n_vec]->SetComponentName(0, "w");
		vectors[n_vec]->SetComponentName(1, "i");
		vectors[n_vec]->SetComponentName(2, "j");
		vectors[n_vec]->SetComponentName(3, "k");
		vectors[n_vec]->SetName("quat");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_QUATpos = n_vec;
		n_vec++;
	}

	if ((fpos[0] >= 0) && (fpos[1] >= 0) && (fpos[2] >= 0)) {
		//vtkErrorMacro(<<"f");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("f");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_fpos = n_vec;
		n_vec++;
	}

	if ((TQpos[0] >= 0) && (TQpos[1] >= 0) && (TQpos[2] >= 0)) {
		//vtkErrorMacro(<<"TQ");
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("tq");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_TQpos = n_vec;
		n_vec++;
	}

	if ((vpos[0] >= 0) && (vpos[1] >= 0) && (vpos[2] >= 0)) {
		//vtkErrorMacro(<<"vel");

		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("v");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_vpos = n_vec;
		n_vec++;
	}

	if ((opos[0] >= 0) && (opos[1] >= 0) && (opos[2] >= 0)) {
		vectors[n_vec] = vtkFloatArray::New();
		vectors[n_vec]->SetNumberOfComponents(3);
		vectors[n_vec]->SetComponentName(0, "X");
		vectors[n_vec]->SetComponentName(1, "Y");
		vectors[n_vec]->SetComponentName(2, "Z");
		vectors[n_vec]->SetName("omega");
		vectors[n_vec]->SetNumberOfTuples(COUNT);
		v_opos = n_vec;
		n_vec++;
	}

	if (RADpos >= 0) {
		//vtkErrorMacro(<<"rad");
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("radius");
		scalars[n_scal]->SetComponentName(0, "radius");
		s_RADpos = n_scal;
		n_scal++;
	}

	if (NUMBONDpos >= 0) {
		scalars[n_scal] = vtkFloatArray::New();
		scalars[n_scal]->SetNumberOfComponents(1);
		scalars[n_scal]->Reset();
		scalars[n_scal]->SetNumberOfTuples(COUNT);
		scalars[n_scal]->SetName("numbonds");
		scalars[n_scal]->SetComponentName(0, "count");
		s_NUMBONDpos = n_scal;
		n_scal++;
	}

	for (int i = 0; i<noo; i++) {
		others[i] = vtkFloatArray::New();
		others[i]->SetNumberOfComponents(1);
		others[i]->Reset();
		others[i]->SetNumberOfTuples(COUNT);
		others[i]->SetName(OTHERname[i].c_str());
	}

	//vtkErrorMacro(<<"got "<<n_scal<<"Scalars");
	//vtkErrorMacro(<<"got "<<n_vec<<"Vectors");

#if defined(WIN32) && defined(MEASURE_TIME)
	QueryPerformanceCounter(&stop);
	MSG(<< "Test3: " << ((double)(stop.QuadPart - start.QuadPart)) / freq);
	start = stop;
#endif

	float val[3];
	val[0] = val[1] = val[2] = 0;
	float vel[3];
	vel[0] = vel[1] = vel[2] = 0;
	float omega[3];
	omega[0] = omega[1] = omega[2] = 0;
	float forces[3];
	forces[0] = forces[1] = forces[2] = 0;
	float ia[3];
	ia[0] = ia[1] = ia[2] = 0;

	float xs[3];
	float xu[3];
	float mu[3];
	float angmom[3];
	double quat[4];
	float tq[3];

	int ic = 0;
	int lc = 0; // Linecounter
	int pc = 0;
	while (this->File->getline(line, sizeof(line))) { //every line
		ic = 0;
		pc = 0;

		char* begin = line;
		char* end;
		float tmp = strtod(begin, &end);
		while (end != begin)
		{
			if (ic == xpos) val[0] = tmp;
			else if (ic == ypos) val[1] = tmp;
			else if (ic == zpos) val[2] = tmp;
			else if (ic == IDpos)		scalars[s_IDpos]->InsertTuple1(lc, (int)tmp);
			else if (ic == TYPEpos)	scalars[s_TYPEpos]->InsertTuple1(lc, (int)tmp);
			else if (ic == RADpos)	scalars[s_RADpos]->InsertTuple1(lc, tmp);
			else if (ic == MOLpos)	scalars[s_MOLpos]->InsertTuple1(lc, tmp);
			else if (ic == MASSpos)	scalars[s_MASSpos]->InsertTuple1(lc, tmp);
			else if (ic == Qpos)		scalars[s_Qpos]->InsertTuple1(lc, tmp);
			else if (ic == DENSpos)	scalars[s_DENSpos]->InsertTuple1(lc, tmp);
			else if (ic == NUMBONDpos) scalars[s_NUMBONDpos]->InsertTuple1(lc, tmp);
			else if (ic == ipos[0]) ia[0] = tmp;
			else if (ic == ipos[1]) ia[1] = tmp;
			else if (ic == ipos[2]) ia[2] = tmp;

			else if (ic == vpos[0]) vel[0] = tmp;
			else if (ic == vpos[1]) vel[1] = tmp;
			else if (ic == vpos[2]) vel[2] = tmp;

			else if (ic == fpos[0]) forces[0] = tmp;
			else if (ic == fpos[1]) forces[1] = tmp;
			else if (ic == fpos[2]) forces[2] = tmp;

			else if (ic == opos[0]) omega[0] = tmp;
			else if (ic == opos[1]) omega[1] = tmp;
			else if (ic == opos[2]) omega[2] = tmp;

			else if (ic == XSpos[0]) xs[0] = tmp;
			else if (ic == XSpos[1]) xs[1] = tmp;
			else if (ic == XSpos[2]) xs[2] = tmp;

			else if (ic == XUpos[0]) xu[0] = tmp;
			else if (ic == XUpos[1]) xu[1] = tmp;
			else if (ic == XUpos[2]) xu[2] = tmp;

			else if (ic == MUpos[0]) mu[0] = tmp;
			else if (ic == MUpos[1]) mu[1] = tmp;
			else if (ic == MUpos[2]) mu[2] = tmp;

			else if (ic == ANGMOMpos[0]) angmom[0] = tmp;
			else if (ic == ANGMOMpos[1]) angmom[1] = tmp;
			else if (ic == ANGMOMpos[2]) angmom[2] = tmp;

			else if (ic == QUATpos[0]) quat[0] = tmp;
			else if (ic == QUATpos[1]) quat[1] = tmp;
			else if (ic == QUATpos[2]) quat[2] = tmp;
			else if (ic == QUATpos[3]) quat[3] = tmp;

			else if (ic == TQpos[0]) tq[0] = tmp;
			else if (ic == TQpos[1]) tq[1] = tmp;
			else if (ic == TQpos[2]) tq[2] = tmp;

			else
			{ //nothing known
				for (int i = 0; i<noo; i++)
					if (ic == OTHERpos[i]) others[i]->InsertTuple1(lc, tmp);
			}

			ic++;
			begin = end;
			tmp = strtod(begin, &end);
		} //item
		points->InsertPoint(lc, val[0], val[1], val[2]);
		if (v_ipos >= 0) vectors[v_ipos]->InsertTuple(lc, ia);
		if (v_vpos >= 0) vectors[v_vpos]->InsertTuple(lc, vel);
		if (v_fpos >= 0) vectors[v_fpos]->InsertTuple(lc, forces);
		if (v_opos >= 0) vectors[v_opos]->InsertTuple(lc, omega);
		if (v_XSpos >= 0) vectors[v_XSpos]->InsertTuple(lc, xs);
		if (v_XUpos >= 0) vectors[v_XUpos]->InsertTuple(lc, xu);
		if (v_MUpos >= 0) vectors[v_MUpos]->InsertTuple(lc, mu);
		if (v_ANGMOMpos >= 0) vectors[v_ANGMOMpos]->InsertTuple(lc, angmom);
		if (v_QUATpos >= 0) vectors[v_QUATpos]->InsertTuple4(lc, quat[0], quat[1], quat[2], quat[3]);
		if (v_TQpos >= 0) vectors[v_TQpos]->InsertTuple(lc, tq);

		lc++;

#if defined(WIN32) && defined(MEASURE_TIME)
		if ((lc % 10000) == 0)
		{
			QueryPerformanceCounter(&stop);
			MSG(<< TS << " - Test3X: " << ((double)(stop.QuadPart - start.QuadPart)) / freq);
			start = stop;
		}
#endif

	}//line
	 //printf("reading done\n");

#if defined(WIN32) && defined(MEASURE_TIME)
	QueryPerformanceCounter(&stop);
	MSG(<< TS << " - Test4: " << ((double)(stop.QuadPart - start.QuadPart)) / freq);
	start = stop;
#endif

	vertices = vtkCellArray::New();
	vertices->Reset();

	this->NumberOfPoints = points->GetNumberOfPoints();
	for (vtkIdType j = 0; j < (vtkIdType)this->NumberOfPoints; ++j)
	{
		vertices->InsertNextCell(1);
		vertices->InsertCellPoint(j);
	}

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the ouptut
	vtkPolyData *myoutput = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	myoutput->SetPoints(points);
	myoutput->SetVerts(vertices);

	if (s_IDpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_IDpos]);
	if (s_TYPEpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_TYPEpos]);
	if (s_RADpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_RADpos]);
	if (s_MASSpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_MASSpos]);
	if (s_MOLpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_MOLpos]);
	if (s_DENSpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_DENSpos]);
	if (s_NUMBONDpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_NUMBONDpos]);
	if (s_Qpos >= 0) myoutput->GetPointData()->AddArray(scalars[s_Qpos]);


	if (v_vpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_vpos]);
	if (v_opos >= 0) myoutput->GetPointData()->AddArray(vectors[v_opos]);
	if (v_fpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_fpos]);
	if (v_ipos >= 0) myoutput->GetPointData()->AddArray(vectors[v_ipos]);

	if (v_XSpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_XSpos]);
	if (v_XUpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_XUpos]);

	if (v_MUpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_MUpos]);
	if (v_ANGMOMpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_ANGMOMpos]);
	if (v_QUATpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_QUATpos]);
	if (v_TQpos >= 0) myoutput->GetPointData()->AddArray(vectors[v_TQpos]);

	/*if (noo>0) {
		for (int i = 0; i<noo; i++) {
			myoutput->GetPointData()->AddArray(others[i]);
		}
	}*/

	vtkIntArray *intValue;
	intValue = vtkIntArray::New();
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

	vtkStringArray *strValue;
	strValue = vtkStringArray::New();
	strValue->SetNumberOfComponents(1);
	strValue->SetName("fname");
	std::string tstr = this->FileName;
	strValue->InsertNextValue(tstr.c_str());
	myoutput->GetFieldData()->AddArray(strValue);
	strValue->Delete();

	if (label.length()>0) {
		strValue = vtkStringArray::New();
		strValue->SetNumberOfComponents(1);
		strValue->SetName("Label");
		strValue->InsertNextValue(label.c_str());
		myoutput->GetFieldData()->AddArray(strValue);
		strValue->Delete();
	}

	myoutput->Modified();

	// Free Memory

	for (int i = 0; i<n_scal; i++) {
		scalars[i]->Delete();
	}
	if (scalars != NULL)
	{
		delete[] scalars;
		scalars = NULL;
	}

	for (int i = 0; i<n_vec; i++) {
		vectors[i]->Delete();
	}
	if (vectors != NULL)
	{
		delete[] vectors;
		vectors = NULL;
	}

	//vtkErrorMacro(<<"delete others");
	for (int i = 0; i<noo; i++) {
		others[i]->Delete();
	}
	if (others != NULL)
	{
		delete[] others;
		others = NULL;
	}

	points->Delete();
	vertices->Delete();

#if defined(WIN32) && defined(MEASURE_TIME)
	QueryPerformanceCounter(&stop);
	MSG(<< "Test5: " << ((double)(stop.QuadPart - start.QuadPart)) / freq);
	start = stop;
#endif

	return 1;
}


int liggghts_reader::RequestInformation(
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

void liggghts_reader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "FileName: "
		<< (this->FileName ? this->FileName : "(NULL)") << endl;
}

int liggghts_reader::CanReadFile(const char *fname)
{
	return 1;
}

