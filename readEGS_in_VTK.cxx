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


// helper class to format slice status message
class StatusMessage {
public:
   static std::string Format(int slice, int maxSlice) {
      std::stringstream tmp;
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
    	  std::cout << "sorting pixel intensity data...\n";
    	  vtkAlgorithm* myProducer {_ProcessedImage->GetProducer()};
    	  //not sure on using port number 1..
    	  vtkDataObject* myData  {myProducer->GetOutputDataObject(0)};

    	  vtkImageData* myImg = myData;




//    	  // Create an image data
//    	  vtkSmartPointer<vtkImageData> imageData =
//    	    vtkSmartPointer<vtkImageData>::New();
//
//    	  // Specify the size of the image data
//    	  imageData->SetDimensions(2,3,1);
//    	#if VTK_MAJOR_VERSION <= 5
//    	  imageData->SetNumberOfScalarComponents(1);
//    	  imageData->SetScalarTypeToDouble();
//    	#else
//    	  imageData->AllocateScalars(VTK_DOUBLE,1);
//    	#endif
//
//    	  int* dims = imageData->GetDimensions();
//    	  // int dims[3]; // can't do this
//
//    	  std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
//
//    	  std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
//    	  std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
//
//    	  // Fill every entry of the image data with "2.0"
//    	  for (int z = 0; z < dims[2]; z++)
//    	    {
//    	    for (int y = 0; y < dims[1]; y++)
//    	      {
//    	      for (int x = 0; x < dims[0]; x++)
//    	        {
//    	        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
//    	        pixel[0] = 2.0;
//    	        }
//    	      }
//    	    }
//
//    	  // Retrieve the entries from the image data and print them to the screen
//    	  for (int z = 0; z < dims[2]; z++)
//    	    {
//    	    for (int y = 0; y < dims[1]; y++)
//    	      {
//    	      for (int x = 0; x < dims[0]; x++)
//    	        {
//    	        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
//    	        // do something with v
//    	        std::cout << pixel[0] << " ";
//    	        }
//    	      std::cout << std::endl;
//    	      }
//    	    std::cout << std::endl;
//    	    }





    	  std::cout << "dumping to text file...\n";
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
   if ( argc != 2 )
   {
      std::cout << "Usage: " << argv[0]
      << " EGSPHANT filename" << std::endl;
      return EXIT_FAILURE;
   }

   uint8_t debugLevel=0;
   egsphantff myObj(debugLevel);
   std::string filename = argv[1];
   myObj.ReadFile(filename);

   std::vector<int16_t> myImg =  myObj.material_arr();
   double dx = myObj.xyz_bound_obj().getXbi(1) -myObj.xyz_bound_obj().getXbi(0);
   double dy = myObj.xyz_bound_obj().getYbi(1) -myObj.xyz_bound_obj().getYbi(0);
   double dz = myObj.xyz_bound_obj().getZbi(1) -myObj.xyz_bound_obj().getZbi(0);

   double x0 = myObj.xyz_bound_obj().getXbi(0);
   double y0 = myObj.xyz_bound_obj().getYbi(0);
   double z0 = myObj.xyz_bound_obj().getZbi(0);

   int32_t nx = myObj.xyz_bound_obj().getNx();
   int32_t ny = myObj.xyz_bound_obj().getNy();
   int32_t nz = myObj.xyz_bound_obj().getNz();



	// Convert the c-style image to a vtkImageData
	vtkSmartPointer<vtkImageImport> imageImport =
	vtkSmartPointer<vtkImageImport>::New();
	imageImport->SetDataSpacing(dx, dy, dz);
	imageImport->SetDataOrigin(x0, y0, z0);
	imageImport->SetWholeExtent(0, nx-1, 0, ny-1, 0, nz-1);
	imageImport->SetDataExtentToWholeExtent();
	imageImport->SetDataScalarTypeToShort();
	imageImport->SetNumberOfScalarComponents(1);
	imageImport->SetImportVoidPointer(myImg.data());
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
