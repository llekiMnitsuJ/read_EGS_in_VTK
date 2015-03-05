//
// This example demonstrates how to read a series of dicom images
// and how to scroll with the mousewheel or the up/down keys
// through all slices
//

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageActor.h>
#include <vtkXMLImageDataWriter.h>

#include <vtkImageGradientMagnitude.h>
#include <vtkImageFlip.h>
#include <vtkAlgorithmOutput.h>
#include <vtkDataObject.h>

// some standard vtk headers
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
// headers needed for this example
#include <vtkImageViewer2.h>
#include <vtkDICOMImageReader.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextMapper.h>
// needed to easily convert int to std::string
#include <sstream>

//we will try to import a material mapping and scroll through it...
#include "egsphantff.hpp"

// i,j,k, are zero based.
// returned 1D index is 0 based
//stride fastest in x direction (third dimension)
template<typename T, typename U, typename S>
S ArrayIndicesX(const U& nx,const U& ny,const U& nz,
		const T& i,const T& j,const T& k,const S& retType){
	S result = i+j*nx+k*nx*ny;
	return result;
}


//use this to hold voxel-based masks
struct mask_output{
	std::string name_;
	int32_t nx_;
	int32_t ny_;
	int32_t nz_;
	std::vector<int16_t> imgData_;
};

void write_mask(const mask_output& mask, const std::string& filename){

	std::ofstream ofs;
	ofs.open(filename.c_str());
	if(ofs.is_open())
	{
		ofs << mask.name_ << "\n";
		ofs << mask.nx_ << " " << mask.ny_ << " " << mask.nz_ << "\n";

		int32_t currIndex=0;
		//x varies fastest, z slowest
		for (int32_t k = 0; k < mask.nz_; ++k)
		{
			for(int32_t j=0; j < mask.ny_; ++j)
			{
				for(int32_t i=0; i < mask.nx_; ++i)
				{
					currIndex = ArrayIndicesX(mask.nx_, mask.ny_, mask.nz_, i, j, k, currIndex);

					ofs << " " << mask.imgData_[currIndex];
				}
				ofs << "\n";
			}
			ofs << "\n";
		}


		ofs.close();
	}
	else
	{
		std::cout << "unable to open " << filename << " for writing!\n";
	}
}

mask_output read_mask(const std::string& filename){
	std::ifstream ifs;
	ifs.open(filename.c_str());

	mask_output result;
	if (ifs.is_open())
	{
		std::cout << "reading: " << filename << "\n";

		std::string temp;
		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		result.name_ = temp;
		std::cout << "name: " << result.name_ << "\n";
		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		std::stringstream ss(temp);
		ss >> result.nx_ >> result.ny_ >> result.nz_;
		std::cout << "nx,ny,nz: " << result.nx_ << "," << result.ny_ << "," << result.nz_ << "\n";


		uint32_t size = result.nx_*result.ny_*result.nz_;
		result.imgData_.clear();
		result.imgData_.resize(size);
		for(uint32_t i=0; i < size-1; ++i) ifs >> result.imgData_[i];
	}
	else
	{
		std::cout << "unable to open file " << filename << "\n";
	}

	return result;
}

struct node{
	uint32_t id_;
	double x_;
	double y_;
	double z_;
};


std::vector<node> bounding_box(const std::vector<node>& arr)
		{
	//low and high values for x,y,z
	double ux=-1E10, lx=1E10, uy=-1E10, ly=1E10, uz=-1E10, lz=1E10;
	for (uint32_t i = 0; i < arr.size(); ++i) {
		if (arr[i].x_ > ux) ux = arr[i].x_;
		if (arr[i].y_ > uy) uy = arr[i].y_;
		if (arr[i].z_ > uz) uz = arr[i].z_;
		if (arr[i].x_ < lx) lx = arr[i].x_;
		if (arr[i].y_ < ly) ly = arr[i].y_;
		if (arr[i].z_ < lz) lz = arr[i].z_;
	}

	node lll;
	lll.id_ = 0;
	lll.x_ = lx;
	lll.y_ = ly;
	lll.z_ = lz;

	node uuu;
	uuu.id_=1;
	uuu.x_ = ux;
	uuu.y_ = uy;
	uuu.z_ = uz;

	std::vector<node> result;
	result.push_back(lll);
	result.push_back(uuu);
	return result;
		}


//Class for sorting
struct jmIndexVal{
	uint16_t i_;
	uint16_t j_;
	uint16_t k_;
	short* valArrPtr_;
};
//Comparison for sorting
bool comp_lt_jmIndexVal(const jmIndexVal& left, const jmIndexVal& right){

	return (left.valArrPtr_[0] < right.valArrPtr_[0]);
}

bool comp_gt_jmIndexVal(const jmIndexVal& left, const jmIndexVal& right){

	return (left.valArrPtr_[0] > right.valArrPtr_[0]);
}
//function for writing pixel location (i,j,k) and intensity value to a file
void writeToFile(const std::string& filename,
		const int& nx, const int& ny, const int& nz,
		const double& x0, const double& y0, const double& z0,
		const double& dx, const double& dy, const double& dz,
		const std::vector<jmIndexVal>& arr){

	std::ofstream ofs{};
	ofs.open(filename.c_str());
	if(ofs.is_open()){

		ofs << nx << " " << ny << " " << nz << "\n";
		ofs << x0 << " " << y0 << " " << z0 << "\n";
		ofs << dx << " " << dy << " " << dz << "\n";

		//val i j k x y z
		for (uint32_t i = 0; i < arr.size(); ++i) {
			ofs << arr[i].valArrPtr_[0] << " "
					<< arr[i].i_ << " "
					<< arr[i].j_ << " "
					<< arr[i].k_ << " "
					<< x0 + (arr[i].i_)*dx << " "
					<< y0 + (arr[i].j_)*dy << " "
					<< z0 + (arr[i].k_)*dz << "\n";
		}

		ofs.close();
	}
	else{
		std::cout << "problem opening " + filename + " for writing.\n";
	}
}


// helper class to format slice status message
class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp{};
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
	}
};


// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	static myVtkInteractorStyleImage* New();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
	vtkImageViewer2* _ImageViewer;
	vtkTextMapper* _StatusMapper;
	vtkAlgorithmOutput* _OriginalImage;
	vtkAlgorithmOutput* _ProcessedImage;

	int _Slice;
	int _MinSlice;
	int _MaxSlice;

public:
	void SetImageViewer(vtkImageViewer2* imageViewer) {
		_ImageViewer = imageViewer;
		_MinSlice = imageViewer->GetSliceMin();
		_MaxSlice = imageViewer->GetSliceMax();
		_Slice = _MinSlice;
		cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
	}

	void SetStatusMapper(vtkTextMapper* statusMapper) {
		_StatusMapper = statusMapper;
	}

	void SetOriginalImage(vtkAlgorithmOutput* myImage){
		_OriginalImage = myImage;
		_ProcessedImage = _OriginalImage;
	}

protected:
	void MoveSliceForward() {
		if(_Slice < _MaxSlice) {
			_Slice += 1;
			cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			_StatusMapper->SetInput(msg.c_str());
			_ImageViewer->Render();
		}
	}

	void MoveSliceBackward() {
		if(_Slice > _MinSlice) {
			_Slice -= 1;
			cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			_StatusMapper->SetInput(msg.c_str());
			_ImageViewer->Render();
		}
	}


	virtual void OnKeyDown() {
		std::string key = this->GetInteractor()->GetKeySym();
		if(key.compare("Up") == 0) {
			//cout << "Up arrow key was pressed." << endl;
			MoveSliceForward();
		}
		else if(key.compare("Down") == 0) {
			//cout << "Down arrow key was pressed." << endl;
			MoveSliceBackward();
		}
		else if(key.compare("g") == 0){
			//Look at the magnitude of the gradient of the image.
			vtkSmartPointer<vtkImageGradientMagnitude> gradMagFilter =
					vtkSmartPointer<vtkImageGradientMagnitude>::New();
			gradMagFilter->SetInputConnection(_ProcessedImage);
			gradMagFilter->SetDimensionality(3);
			gradMagFilter->Update();
			_ProcessedImage = gradMagFilter->GetOutputPort();

			_ImageViewer->SetColorWindow(10);
			_ImageViewer->SetColorLevel(0.5*10);
			_ImageViewer->SetInputConnection(_ProcessedImage);
			_ImageViewer->Render();

		}
		else if(key.compare("o") == 0){
			_ProcessedImage = _OriginalImage;
			_ImageViewer->SetInputConnection(_ProcessedImage);
			_ImageViewer->SetColorWindow(10);
			_ImageViewer->SetColorLevel(0.5*10);
			_ImageViewer->Render();
		}
		else if(key.compare("y") == 0){
			//TODO: this doesnt work so well once I translate the image...
			//flip the image since it is upside down on screen
			vtkSmartPointer<vtkImageFlip> flipYFilter =
					vtkSmartPointer<vtkImageFlip>::New();
			flipYFilter->SetFilteredAxis(1); // flip y axis
			flipYFilter->SetOutputDimensionality(3);
			flipYFilter->SetInputConnection(_ProcessedImage);
			flipYFilter->Update();
			_ProcessedImage = flipYFilter->GetOutputPort();

			_ImageViewer->SetColorWindow(10);
			_ImageViewer->SetColorLevel(0.5*10);
			_ImageViewer->SetInputConnection(_ProcessedImage);
			_ImageViewer->Render();
		}
		else if(key.compare("d") == 0){
			//dump pixel values sorted by intensity to a file


			vtkAlgorithm* myProducer {_ProcessedImage->GetProducer()};
			//not sure on using port number 0..
			vtkDataObject* myData  {myProducer->GetOutputDataObject(0)};

			vtkImageData* imageData = vtkImageData::SafeDownCast(myData);


			int* dims = imageData->GetDimensions();
			double* origin = imageData->GetOrigin();
			double* spacing = imageData->GetSpacing();
			std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
			std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
			std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
			std::cout << "Origin: " << "(" << origin[0] << "," << origin[1] << "," << origin[2] << ")\n";
			std::cout << "Spacing: " << "(" << spacing[0] << "," << spacing[1] << "," << spacing[2] << ")\n";
			std::cout << "sorting pixel intensity data...\n";




			std::vector<jmIndexVal> myIndices{};
			// copy every pixel and value to a temporary datastructure
			for (int z = 0; z < dims[2]; z++)
			{
				for (int y = 0; y < dims[1]; y++)
				{
					for (int x = 0; x < dims[0]; x++)
					{
						jmIndexVal temp{};
						temp.i_ = x; temp.j_ = y; temp.k_ = z;
						temp.valArrPtr_ = static_cast<short*>(imageData->GetScalarPointer(x,y,z));
						myIndices.push_back(temp);
					}
				}
			}


			std::sort(myIndices.begin(), myIndices.end(), comp_gt_jmIndexVal);
			std::cout << "dumping to text file...\n";
			writeToFile("sortedPoints.txt",
					dims[0], dims[1],dims[2],
					origin[0], origin[1], origin[2],
					spacing[0], spacing[1], spacing[2],
					myIndices);



		}

		// forward event
		vtkInteractorStyleImage::OnKeyDown();
	}


	virtual void OnMouseWheelForward() {
		//std::cout << "Scrolled mouse wheel forward." << std::endl;
		MoveSliceForward();
		// don't forward events, otherwise the image will be zoomed
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelForward();
	}


	virtual void OnMouseWheelBackward() {
		//std::cout << "Scrolled mouse wheel backward." << std::endl;
		if(_Slice > _MinSlice) {
			MoveSliceBackward();
		}
		// don't forward events, otherwise the image will be zoomed
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelBackward();
	}
};

vtkStandardNewMacro(myVtkInteractorStyleImage);


int main(int argc, char* argv[])
{
	// Verify input arguments
	if ( argc != 3 )
	{
		std::cout << "Usage: " << argv[0]
									   << " EGSPHANT_filename"
									   << " body_mask_file" << std::endl;
		return EXIT_FAILURE;
	}

	uint8_t debugLevel=0;
	std::string filename = argv[1];
	std::string maskFilename = argv[2];

	egsphantff myObj(debugLevel);
	myObj.ReadFile(filename);

	std::vector<double> xb = myObj.xyz_bound_obj().xb_arr();
	std::vector<double> yb = myObj.xyz_bound_obj().yb_arr();
	std::vector<double> zb = myObj.xyz_bound_obj().zb_arr();

	std::vector<int16_t> myImg =  myObj.material_arr();
	double dx {xb[1] - xb[0]};
	double dy {yb[1] - yb[0]};
	double dz {zb[1] - zb[0]};

	double x0 = myObj.xyz_bound_obj().getXbi(0);
	double y0 = myObj.xyz_bound_obj().getYbi(0);
	double z0 = myObj.xyz_bound_obj().getZbi(0);

	int32_t nx = myObj.xyz_bound_obj().getNx();
	int32_t ny = myObj.xyz_bound_obj().getNy();
	int32_t nz = myObj.xyz_bound_obj().getNz();


	std::cout << "reading mask file: " << maskFilename << "\n";
	mask_output body_mask = read_mask(maskFilename); //this has 0s and 1s
	std::vector<node> myNodes{}; //this only has 1s

	uint32_t nodeID=0, currIndex=0;
	node myNode{};
	for (int32_t i = 0; i < body_mask.nx_; ++i)
		for (int32_t j = 0; j < body_mask.ny_; ++j)
			for (int32_t k = 0; k < body_mask.nz_; ++k)
			{
				currIndex = ArrayIndicesX(body_mask.nx_,
						body_mask.ny_,
						body_mask.nz_, i,j,k,currIndex);
				if(body_mask.imgData_[currIndex] != 0)
				{
					myNode.id_ = nodeID;
					myNode.x_ = (xb[i]+xb[i+1])/2.;
					myNode.y_ = (yb[j]+yb[j+1])/2.;
					myNode.z_ = (zb[k]+zb[k+1])/2.;
					myNodes.push_back(myNode);
					nodeID++;
				}
				else
				{
					//zero out data in egsphant material
					myImg[currIndex] = 1; //1 is min value
				}
			}

	std::vector<node> bb = bounding_box(myNodes);
	int startx = floor((bb[0].x_ - xb[0])/dx);
	int starty = floor((bb[0].y_ - yb[0])/dy);
	int startz = floor((bb[0].z_ - zb[0])/dz);
	int endx = floor((bb[1].x_ - xb[0])/dx);
	int endy = floor((bb[1].y_ - yb[0])/dy);
	int endz = floor((bb[1].z_ - zb[0])/dz);

	//just pad a little more...
	if (startx > 0) startx--;
	if (starty > 0) starty--;
	if (startz > 0) startz--;
	if(endx < (nx-1)) endx++;
	if(endy < (ny-1)) endy++;
	if(endz < (nz-1)) endz++;

	BOOST_ASSERT(startx >= 0);
	BOOST_ASSERT(starty >= 0);
	BOOST_ASSERT(startz >= 0);
	BOOST_ASSERT(endx < nx);
	BOOST_ASSERT(endy < ny);
	BOOST_ASSERT(endz < nz);

	std::vector<int16_t> croppedImg{};
	for (int32_t k = startz; k <= endz; ++k)
		for (int32_t j = starty; j <= endy; ++j)
			for (int32_t i = startx; i <= endx; ++i){
				currIndex = ArrayIndicesX(nx,ny,nz,i,j,k,currIndex);
				croppedImg.push_back(myImg[currIndex]);
			}


	std::cout << "cropping to ["<< startx << ":" << endx << ","
			                    << starty << ":" << endy << ","
								<< startz << ":" << endz << "]\n";

	std::cout << "cropped image has " << croppedImg.size() << " elements\n";
	//sub array has new size...
	int newNx {endx - startx + 1};
	int newNy {endy - starty + 1};
	int newNz {endz - startz + 1};
	std::cout << "cropped nx,ny,nz: " << newNx << " " << newNy << " " << newNz << "\n";

	double newX0 {(xb[startx] +xb[startx+1])/2};
	double newY0 {(yb[starty] +yb[starty+1])/2};
	double newZ0 {(zb[startz] +zb[startz+1])/2};
	std::cout << "cropped x0,y0,z0: " << newX0 << " " << newY0 << " " << newZ0 << "\n";

	// Convert the c-style image to a vtkImageData
	vtkSmartPointer<vtkImageImport> imageImport =
			vtkSmartPointer<vtkImageImport>::New();
	imageImport->SetDataSpacing(dx, dy, dz);
	imageImport->SetDataOrigin(newX0,newY0,newZ0);
	imageImport->SetWholeExtent(0, newNx-1, 0, newNy-1, 0, newNz-1);
	imageImport->SetDataExtentToWholeExtent();
	imageImport->SetDataScalarTypeToShort();
	imageImport->SetNumberOfScalarComponents(1);
	imageImport->SetImportVoidPointer(croppedImg.data());
	imageImport->Update();










		   // Visualize
		   vtkSmartPointer<vtkImageViewer2> imageViewer =
				   vtkSmartPointer<vtkImageViewer2>::New();
//imageViewer->SetInputConnection(imageImport->GetOutputPort());
imageViewer->SetInputConnection(imageImport->GetOutputPort());


// slice status message
vtkSmartPointer<vtkTextProperty> sliceTextProp = vtkSmartPointer<vtkTextProperty>::New();
sliceTextProp->SetFontFamilyToCourier();
sliceTextProp->SetFontSize(10);
sliceTextProp->SetVerticalJustificationToBottom();
sliceTextProp->SetJustificationToLeft();

vtkSmartPointer<vtkTextMapper> sliceTextMapper = vtkSmartPointer<vtkTextMapper>::New();
std::string msg = StatusMessage::Format(imageViewer->GetSliceMin(), imageViewer->GetSliceMax());
sliceTextMapper->SetInput(msg.c_str());
sliceTextMapper->SetTextProperty(sliceTextProp);

vtkSmartPointer<vtkActor2D> sliceTextActor = vtkSmartPointer<vtkActor2D>::New();
sliceTextActor->SetMapper(sliceTextMapper);
sliceTextActor->SetPosition(15, 10);

// usage hint message
vtkSmartPointer<vtkTextProperty> usageTextProp = vtkSmartPointer<vtkTextProperty>::New();
usageTextProp->SetFontFamilyToCourier();
usageTextProp->SetFontSize(8);
usageTextProp->SetVerticalJustificationToTop();
usageTextProp->SetJustificationToLeft();

vtkSmartPointer<vtkTextMapper> usageTextMapper = vtkSmartPointer<vtkTextMapper>::New();
usageTextMapper->SetInput("- Slice with mouse wheel\n  or Up/Down-Key\n- Zoom with pressed right\n  mouse button while dragging");
usageTextMapper->SetTextProperty(usageTextProp);

vtkSmartPointer<vtkActor2D> usageTextActor = vtkSmartPointer<vtkActor2D>::New();
usageTextActor->SetMapper(usageTextMapper);
usageTextActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
usageTextActor->GetPositionCoordinate()->SetValue( 0.05, 0.95);

// create an interactor with our own style (inherit from vtkInteractorStyleImage)
// in order to catch mousewheel and key events
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle =
		vtkSmartPointer<myVtkInteractorStyleImage>::New();

// make imageviewer2 and sliceTextMapper visible to our interactorstyle
// to enable slice status message updates when scrolling through the slices
myInteractorStyle->SetImageViewer(imageViewer);
myInteractorStyle->SetStatusMapper(sliceTextMapper);
myInteractorStyle->SetOriginalImage(imageImport->GetOutputPort());



imageViewer->SetupInteractor(renderWindowInteractor);
// make the interactor use our own interactorstyle
// cause SetupInteractor() is defining it's own default interatorstyle
// this must be called after SetupInteractor()
renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
// add slice status message and usage hint message to the renderer
imageViewer->GetRenderer()->AddActor2D(sliceTextActor);
imageViewer->GetRenderer()->AddActor2D(usageTextActor);

// initialize rendering and interaction
//imageViewer->GetRenderWindow()->SetSize(400, 300);
//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
imageViewer->Render();
imageViewer->GetRenderer()->ResetCamera();
//set the window level appropriate for material indices
imageViewer->SetColorWindow(10);
imageViewer->SetColorLevel(0.5*10);
imageViewer->Render();

renderWindowInteractor->Start();

return EXIT_SUCCESS;
}
