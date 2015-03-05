/*
 * nodeselection.h
 *
 *  Created on: Mar 4, 2015
 *      Author: jkmikell
 */

#ifndef NODE_SELECTION_H_
#define NODE_SELECTION_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <string>
#include <boost/assert.hpp>

template <typename Z>
struct nodeIndexVal{
	int i_;
	int j_;
	int k_;
	Z val_;
};

template <typename Z>
inline bool comp_nodeIndexVal(const nodeIndexVal<Z> & left, const nodeIndexVal<Z>& right)
{
	return (left.val_ < right.val_);
}

template <typename Z>
class node_selection {

public:
	node_selection();
	virtual ~node_selection();
	void ReadFile(const std::string& filename);
	void WriteFile(const std::string& filename);
	node_selection<double> getProbability();
	node_selection<double> CreateListOfVertices(
			const double& radius, const uint32_t maxPoints=0);

	// i,j,k, are zero based.
	// returned 1D index is 0 based
	//stride fastest in x direction (third dimension)
	template<typename T, typename U, typename S>
	S ArrayIndicesX(const U& nx,const U& ny,const U& nz,
			const T& i,const T& j,const T& k,const S& retType) const {
		S result = i+j*nx+k*nx*ny;
		return result;
	}

	template<typename T, typename U, typename S>
	S GetXLoc(const U& nx, const U& ny, const U& nz,
			const T& i, const T& j, const T& k, const S& retType) const {
		S result = x0_ + i*(nx-1)*dx_;
		return result;
	}

	template<typename T, typename U, typename S>
	S GetYLoc(const U& nx, const U& ny, const U& nz,
			const T& i, const T& j, const T& k, const S& retType) const {
		S result = y0_ + j*(ny-1)*dy_;
		return result;
	}

	template<typename T, typename U, typename S>
	S GetZLoc(const U& nx, const U& ny, const U& nz,
			const T& i, const T& j, const T& k, const S& retType) const {
		S result = z0_ + k*(nz-1)*dz_;
		return result;
	}

	int nx_;
	int ny_;
	int nz_;
	double x0_;
	double y0_;
	double z0_;
	double dx_;
	double dy_;
	double dz_;
	std::vector< nodeIndexVal<Z> > imgArr_;
};

template<typename Z>
inline void node_selection<Z>::ReadFile(const std::string& filename) {
	std::ifstream ifs{};
	ifs.open(filename.c_str());
	if (ifs.is_open()){
		std::cout << "reading: " << filename << "\n";

		std::string temp{};
		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		std::stringstream ss{};
		ss.clear();
		ss.str(temp);
		ss >> nx_ >> ny_ >> nz_;
		std::cout << "nx,ny,nz: " << nx_ << "," << ny_ << "," << nz_  << "\n";

		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		ss.clear();
		ss.str(temp);
		ss >> x0_ >> y0_ >> z0_;
		std::cout << "x0,y0,z0: " << x0_ << "," << y0_ << "," << z0_  << "\n";

		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		ss.clear();
		ss.str(temp);
		ss >> dx_ >> dy_ >> dz_;
		std::cout << "dx,dy,dz: " << dx_ << "," << dy_ << "," << dz_  << "\n";


		uint32_t size = nx_*ny_*nz_;
		imgArr_.clear();
		imgArr_.resize(size);
		double tempdoub{};
		for(uint32_t i=0; i < size; ++i){
			if (!std::getline(ifs, temp)) std::cout << "problem reading " <<
					filename << " at i= " << i << "\n";
			ss.clear();
			ss.str(temp);
			ss >> imgArr_[i].val_ >>
			imgArr_[i].i_ >> imgArr_[i].j_ >> imgArr_[i].k_ >>
			tempdoub >> tempdoub >> tempdoub;
		}
		ifs.close();
	}
	else{
		std::cout << "unable to open file for reading " << filename << "\n";
	}

}

template<typename Z>
inline void node_selection<Z>::WriteFile(const std::string& filename) {

	std::ofstream ofs{};
	ofs.open(filename.c_str());
	if(ofs.is_open()){

		std::cout << "writing " << filename << "\n";
		ofs << nx_ << " " << ny_ << " " << nz_ << "\n";
		ofs << x0_ << " " << y0_ << " " << z0_ << "\n";
		ofs << dx_ << " " << dy_ << " " << dz_ << "\n";

		//val i j k x y z
		for (uint32_t i = 0; i < imgArr_.size(); ++i) {
			ofs << imgArr_[i].val_ << " "
					<< imgArr_[i].i_ << " "
					<< imgArr_[i].j_ << " "
					<< imgArr_[i].k_ << " "
					<< x0_ + imgArr_[i].i_*dx_ << " "
					<< y0_ + imgArr_[i].j_*dy_ << " "
					<< z0_ + imgArr_[i].k_*dz_ << "\n";
		}

		ofs.close();
	}
	else{
		std::cout << "problem opening " + filename + " for writing.\n";
	}

}

template<typename Z>
inline node_selection<Z>::node_selection():nx_(1), ny_(1), nz_(1),
x0_(0), y0_(0), z0_(0),
dx_(1), dy_(1), dz_(1){

	imgArr_ = std::vector< nodeIndexVal<Z> >{};
	imgArr_.clear();
	imgArr_.resize(1);
}

template<typename Z>
inline node_selection<Z>::~node_selection() {
}

template<typename Z>
inline node_selection<double> node_selection<Z>::getProbability() {
	//here we just normalize the image data so the sum of all voxels is 1.
	double total{0}, checkTotal{0};
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		total += imgArr_[i].val_;
	}
	std::cout << "elements in original image: " << imgArr_.size() << "\n";
	std::cout << "sum of original image: " << total << "\n";

	auto myMin = std::min_element(imgArr_.begin(), imgArr_.end(), comp_nodeIndexVal<Z>);
	auto myMax = std::max_element(imgArr_.begin(), imgArr_.end(), comp_nodeIndexVal<Z>);

	std::cout << "original values range from " << (*myMin).val_ << " to " << (*myMax).val_ << "\n";

	node_selection<double> result{};

	result.nx_ = nx_;
	result.ny_ = ny_;
	result.nz_ = nz_;
	result.x0_ = x0_;
	result.y0_ = y0_;
	result.z0_ = z0_;
	result.dx_ = dx_;
	result.dy_ = dy_;
	result.dz_ = dz_;
	result.imgArr_.clear();
	result.imgArr_.resize(imgArr_.size());

	BOOST_ASSERT(result.imgArr_.size() == imgArr_.size());
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		result.imgArr_[i].i_ = imgArr_[i].i_;
		result.imgArr_[i].j_ = imgArr_[i].j_;
		result.imgArr_[i].k_ = imgArr_[i].k_;
		result.imgArr_[i].val_ =  imgArr_[i].val_/total;
		checkTotal+=result.imgArr_[i].val_;
	}


	std::cout << "elements in normalized image: " << result.imgArr_.size() << "\n";
	std::cout << "total of normalized image: " << checkTotal << "\n";

	return result;

}


//this assumes the current imgArr is a sorted probability node_selection
// see readEGS_in_VTK.cxx
//
//e.g. from
/*
 * int main(int argc, char* argv[])
{
	// Verify input arguments
	if ( argc != 3 )
	{
		std::cout << "Usage: " << argv[0]
									   << " sorted_voxel_file"
									   << " output_probability_filename" << std::endl;
		return 1;
	}
	uint8_t debugLevel=0;
	std::string sortedFilename = argv[1];
	std::string probFilename = argv[2];
	node_selection<int16_t> myNodeObj{};
	myNodeObj.ReadFile(sortedFilename);
	std::cout << "myNodeObj imgArr size: " << myNodeObj.imgArr_.size() << "\n";
	node_selection<double> myProbObj = myNodeObj.getProbability(); ********<-----myProbObj is now ready for this function

	return 0;
}
 *
 */
//
template<typename Z>
inline node_selection<double> node_selection<Z>::CreateListOfVertices(
		const double& radius,
		const uint32_t maxPoints){

	uint32_t myMaxPoints{maxPoints};
	if  (myMaxPoints == 0){
		myMaxPoints = imgArr_.size();
	}
	std::cout << "maximum number of points to generate: " << myMaxPoints << "\n";

	uint32_t nonZeroProbPts{0};

	//generate an image to hold the Pijk
	std::vector<double> pijk{};
	pijk.clear();
	pijk.resize(imgArr_.size());
	uint32_t currIndex=0;
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		currIndex= ArrayIndicesX(nx_, ny_, nz_,
				imgArr_[i].i_,imgArr_[i].j_,imgArr_[i].k_, currIndex);
		pijk[currIndex] = imgArr_[i].val_;
		if (imgArr_[i].val_ != 0.) nonZeroProbPts++;
	}

	std::cout << " fraction of nonzero prob pts: " << nonZeroProbPts << "/" << imgArr_.size() << "\n";
	std::cout << " throwing away zero prob pts...\n";
	imgArr_.resize(nonZeroProbPts);
	nonZeroProbPts=0;
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
			if (imgArr_[i].val_ != 0.) nonZeroProbPts++;
		}
	std::cout << " fraction of nonzero prob pts: " << nonZeroProbPts << "/" << imgArr_.size() << "\n";

	node_selection<double> result{};

	//add max probability to result



	return result;
}

#endif /* NODE_SELECTION_H_ */
