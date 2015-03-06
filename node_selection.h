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
inline bool comp_lt_nodeIndexVal(const nodeIndexVal<Z>& left, const nodeIndexVal<Z>& right)
{
	return (left.val_ < right.val_);
}

template <typename Z>
inline bool comp_gt_nodeIndexVal(const nodeIndexVal<Z>& left, const nodeIndexVal<Z>& right)
{
	return (left.val_ > right.val_);
}


template <typename Z>
std::ostream& operator<<(std::ostream& os, const nodeIndexVal<Z>& o){
	os << o.val_ << " "  << o.i_ << " " << o.j_ << " " << o.k_;
	return(os);
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (std::size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});

  return idx;
}
template <typename T>
double calcRadius(const T& i0, const T& j0, const T& k0,
		const T& i, const T&j, const T& k,
		const double& dx, const double& dy, const double& dz)
{
	return sqrt(pow(dx*(i-i0), 2) + pow(dy*(j-j0),2) + pow(dz*(k-k0),2));
}

template <typename Z>
double calcRadius(const nodeIndexVal<Z>& o, const nodeIndexVal<Z>& p,
		const double& dx, const double& dy, const double& dz)
{
	return calcRadius(o.i_,o.j_,o.k_,
			p.i_,p.j_,p.k_,
			dx,dy,dz);
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
	void updateProbabilities(const node_selection<double>& ineligibleNodes,
			const double& radius);
	double updateProbabilities(const nodeIndexVal<double>& centerNode,
			const double& radius, std::vector<double>& probArr);
	void updateAndSort(std::vector<double>& probArr, std::vector<bool>& okToAddNode);


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

		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		ss.clear();
		ss.str(temp);
		uint32_t imgArrSize{0};
		ss >> imgArrSize;
		std::cout << "elements in imgArr: " << imgArrSize << "\n";

		imgArr_.clear();
		imgArr_.resize(imgArrSize);
		double tempdoub{};
		for(uint32_t i=0; i < imgArrSize; ++i){
			if (!std::getline(ifs, temp)) std::cout << "problem reading " <<
					filename << " at i= " << i << "\n";
			ss.clear();
			ss.str(temp);
			ss >> imgArr_[i].val_ >>
			imgArr_[i].i_ >> imgArr_[i].j_ >> imgArr_[i].k_ >>
			tempdoub >> tempdoub >> tempdoub;
		}
		ifs.close();
		std::cout << "finsihed reading file " << filename << "\n";
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
		ofs << imgArr_.size() << "\n";

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

	auto myMin = std::min_element(imgArr_.begin(), imgArr_.end(), comp_lt_nodeIndexVal<Z>);
	auto myMax = std::max_element(imgArr_.begin(), imgArr_.end(), comp_lt_nodeIndexVal<Z>);

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
node_selection<double> node_selection<Z>::CreateListOfVertices(
		const double& radius,
		const uint32_t maxPoints){

	uint32_t myMaxPoints{maxPoints};
	if  (myMaxPoints == 0){
		myMaxPoints = imgArr_.size();
	}
	std::cout << "maximum number of points to generate: " << myMaxPoints << "\n";

	uint32_t nonZeroProbPts{0};
	uint32_t currIndex{0};
	for (uint32_t i = 0; i < imgArr_.size(); ++i) if (imgArr_[i].val_ != 0.) nonZeroProbPts++;
	std::cout << " fraction of nonzero prob pts: " << nonZeroProbPts << "/" << imgArr_.size() << "\n";
	std::cout << " throwing away zero prob pts...\n";
	imgArr_.resize(nonZeroProbPts);
	nonZeroProbPts=0;
	for (uint32_t i = 0; i < imgArr_.size(); ++i) if (imgArr_[i].val_ != 0.) nonZeroProbPts++;
	std::cout << " fraction of nonzero prob pts: " << nonZeroProbPts << "/" << imgArr_.size() << "\n";

	node_selection<double> result{};
	result.x0_ = x0_;
	result.y0_ = y0_;
	result.z0_ = z0_;
	result.dx_ = dx_;
	result.dy_ = dy_;
	result.dz_ = dz_;
	result.nx_ = nx_;
	result.ny_ = ny_;
	result.nz_ = nz_;

	//initialize array
	result.imgArr_.clear();
	result.imgArr_.resize(imgArr_.size());
	for (uint32_t i = 0; i < result.imgArr_.size(); ++i) {
		result.imgArr_[i].i_ = -1;
		result.imgArr_[i].j_ = -1;
		result.imgArr_[i].k_ = -1;
		result.imgArr_[i].val_ = -1;
	}

	//fill a dummy probArr so we can move easily through space
	std::vector<double> probArr(nx_*ny_*nz_, 0.0);
	//now update the dummyArr with current true probs
	currIndex=0;
	double totalProb{0};
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		currIndex = ArrayIndicesX(nx_, ny_, nz_,
				imgArr_[i].i_,imgArr_[i].j_,imgArr_[i].k_, currIndex);
		probArr[currIndex] = imgArr_[i].val_;
		totalProb+=probArr[currIndex];
	}
	std::cout << "total prob before updates...: " << totalProb << "\n";

	//we will keep the imgArr_ sorted, just stop one short to avoid bounds error
	double nextVal{0}; //this will hold the probability for the next entry
	double maxProbUpdate{0}, currentMax{0};
	uint32_t dummyIndex{0};
	std::vector<bool> okToAddNode(imgArr_.size(), true);
	for (uint32_t i = 0; i < imgArr_.size()-1; ++i) {
		if(okToAddNode[i]){
			std::cout << "updating probability around " << imgArr_[i] << "\n";
			nextVal = imgArr_[i+1].val_;
			//update probability and track maximum returned probability
			maxProbUpdate = updateProbabilities(imgArr_[i], radius, probArr);
			//keep track of maximum updated value
			if (maxProbUpdate > currentMax) currentMax = maxProbUpdate;
			if (currentMax > nextVal){
				result.imgArr_[dummyIndex] = imgArr_[i];
				okToAddNode[i] = false;
				i=0;
				updateAndSort(probArr, okToAddNode);
			}
			else{
				result.imgArr_[dummyIndex] = imgArr_[i];
				okToAddNode[i] = false;
			}
		}
	}



	/*
	//create a type to hold state of node
	// eligible (it can still be added to result)
	// ineligible (we have visited and added vertex to result, so it can no longer be added)
	// temp_ineligible (we have visited and marked it too close to another node, but it can become eligible on a second pass)
	enum jBool:uint8_t {ELIGIBLE=0, INELIGIBLE=1, TEMP_INELIGIBLE=2};
	std::vector<jBool> eligibleNode(imgArr_.size(), jBool::ELIGIBLE);
	std::cout << "eligible node size: " << eligibleNode.size() << "\n";

	//the node at top of list automatically goes in since it is sorted...
	std::cout << "adding first node: " << imgArr_[0] << "\n";
	uint32_t dummyIndex=0;
	result.imgArr_[dummyIndex] = imgArr_[0];
	dummyIndex++;
	eligibleNode[0] = jBool::INELIGIBLE;

	//update the probabilities of its neighbors


	nodeIndexVal<double> nextNode{};
	//go through list and select points that are > radius away from
	for (uint32_t i = 1; i < imgArr_.size(); ++i) {
		std::cout << "\ntesting " << i << "/" << imgArr_.size() << " to see if we should add to result: " << imgArr_[i] << "\n";
		if(eligibleNode[i]==jBool::ELIGIBLE)
		{
			std::cout << "initially eligible\n";
			nextNode = imgArr_[i];
			//check if nextNode is within radius of any result nodes and mark invalid if it is...
			double myCalcRadius=0;
			for(uint32_t j=0; j < result.imgArr_.size(); ++j){
				myCalcRadius = calcRadius(nextNode, result.imgArr_[j], dx_, dy_, dz_);
				if( myCalcRadius <= radius){
					eligibleNode[i] = jBool::TEMP_INELIGIBLE;
					std::cout << "marked temporary ineligible: within radius (" << myCalcRadius << " < " << radius << ") of a result node\n";
					break;
				}
			}
			//if still eligible then add to result and mark ineligible
			if(eligibleNode[i] == jBool::ELIGIBLE){
				std::cout << "adding node to result: " << nextNode << "\n";
				result.imgArr_[dummyIndex] = nextNode;
				dummyIndex++;
				eligibleNode[i] = jBool::INELIGIBLE;
				std::cout << "total result nodes: " << dummyIndex << "\n";
			}
		}
		else{
			std::cout << "was ineligible\n";
		}
	}

	//now trim result index
	result.imgArr_.resize(dummyIndex);

	//update the internal imgArr_ values by moving ELIGIBLE and TEMP_INELIGIBLE to ELIGIBLE.
	uint32_t nE{0}, nIE{0}, ntempIE{0};
	std::vector<uint32_t> ntempIE_indices{};
	for (uint32_t i=0; i < eligibleNode.size(); ++i){
		switch (eligibleNode[i]) {
		case jBool::ELIGIBLE:
			nE++;
			ntempIE_indices.push_back(i);
			break;
		case jBool::INELIGIBLE:
			nIE++;
			break;
		case jBool::TEMP_INELIGIBLE:
			ntempIE++;
			ntempIE_indices.push_back(i);
			break;
		default:
			break;
		}
	}

	std::cout << "after pass there are " << nE << " eligible, "
			<< nIE << " ineligible, and " << ntempIE << " temp ineligible\n";

	std::cout << "updating internal imgArr_ with temp_ineligible members...\n";
	for (uint32_t i = 0; i < ntempIE_indices.size(); ++i) {
		imgArr_[i] = imgArr_[ntempIE_indices[i]];
	}


	std::cout << "trimming internal imgArr_ to " << ntempIE_indices.size() << "\n";
	imgArr_.resize(ntempIE_indices.size());

	std::cout << "updating probabilities...\n";
	updateProbabilities(result, radius);
	std::sort(imgArr_.begin(), imgArr_.end(), comp_gt_nodeIndexVal<double>);
	 */
	return result;
}

//this function iterates through the ineligbleNodes and applies
// a spatial function to adjust neighboring probabilities.
// this believes that the current object imgArr_ only contains eligible or temp_ineligible nodes.
template<typename Z>
inline void node_selection<Z>::updateProbabilities(
		const node_selection<double>& ineligibleNodes, const double& radius) {



	//fill a dummy probArr so we can move easily through space
	std::vector<double> dummyArr(nx_*ny_*nz_, 0.0);

	//now update the dummyArr with current true probs
	uint32_t currIndex=0;
	double totalProb{0};
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		currIndex = ArrayIndicesX(nx_, ny_, nz_,
				imgArr_[i].i_,imgArr_[i].j_,imgArr_[i].k_, currIndex);
		dummyArr[currIndex] = imgArr_[i].val_;
		totalProb+=dummyArr[currIndex];
	}
	std::cout << "total prob before updates...: " << totalProb << "\n";

	//now go through the ineligibleNodes and decrease probabilities of their neighbors
	int gMinI, gMaxI, gMinJ, gMaxJ, gMinK, gMaxK; //global min/max indices
	int ti, tj, tk; //temp indices
	//find bounding integers for extent of radius
	int32_t const deltaI {ceil(radius/dx_)};
	int32_t const deltaJ {ceil(radius/dy_)};
	int32_t const deltaK {ceil(radius/dz_)};
	for (uint32_t i = 0; i < ineligibleNodes.imgArr_.size(); ++i) {

		ti = ineligibleNodes.imgArr_[i].i_;
		tj = ineligibleNodes.imgArr_[i].j_;
		tk = ineligibleNodes.imgArr_[i].k_;


		gMinI = ti - deltaI;
		gMinI = (gMinI < 0) ? 0 : gMinI;
		gMinJ = tj - deltaJ;
		gMinJ = (gMinJ < 0) ? 0 : gMinJ;
		gMinK = tk - deltaK;
		gMinK = (gMinK < 0) ? 0 : gMinK;

		gMaxI = ti + deltaI;
		gMaxI = (gMaxI >= nx_) ? nx_ - 1 : gMaxI;
		gMaxJ = tj + deltaJ;
		gMaxJ = (gMaxJ >= ny_) ? ny_ - 1 : gMaxJ;
		gMaxK = tk + deltaK;
		gMaxK = (gMaxK >= nz_) ? nz_ - 1 : gMaxK;

		double myRadius=0;
		for (int ii = gMinI; ii <= gMaxI; ++ii)
			for(int jj = gMinJ; jj <= gMaxJ; ++jj)
				for(int kk = gMinK; kk <= gMaxK; ++kk){
					//check that its nonzero
					currIndex = ArrayIndicesX(nx_,ny_, nz_,
							ii,jj,kk, currIndex);
					if (dummyArr[currIndex] > 0) {
						myRadius = calcRadius(ti,tj,tk, ii,jj,kk, dx_, dy_,dz_);
						if (myRadius < radius){
							dummyArr[currIndex] *= pow(myRadius/radius, 2); //just using a quadratic, by definition will vary between 0 an and 1
						}

					}

				}
	}


	totalProb=0;
	for(uint32_t i=0; i < dummyArr.size(); ++i) totalProb+=dummyArr[i];
	std::cout << "total prob after updates...: " << totalProb << "\n";

	//now update the dummyArr with current true probs
	std::cout << "updating the nodes imgArr_ ...\n";
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		currIndex = ArrayIndicesX(nx_, ny_, nz_,
				imgArr_[i].i_,imgArr_[i].j_,imgArr_[i].k_, currIndex);
		imgArr_[i].val_ = dummyArr[currIndex];
	}

	std::cout << "sorting nodes by prob...\n";
	std::sort(imgArr_.begin(), imgArr_.end(), comp_gt_nodeIndexVal<double>);


}

//use this to update a generated probArr around the centerNode
//it decreases probabilities of centerNodes neighbors
//it returns the highest modified probability which can be useful for determining when to sort the original array.
template<typename Z>
inline double node_selection<Z>::updateProbabilities(
		const nodeIndexVal<double>& centerNode, const double& radius,
		std::vector<double>& probArr) {

	//find bounding integers for extent of radius
	const int deltaI {ceil(radius/dx_)};
	const int deltaJ {ceil(radius/dy_)};
	const int deltaK {ceil(radius/dz_)};

	const int ti {centerNode.i_};
	const int tj {centerNode.j_};
	const int tk {centerNode.k_};


	//global min/max indices
	int gMinI {ti - deltaI};
	gMinI = (gMinI < 0) ? 0 : gMinI;
	int gMinJ {tj - deltaJ};
	gMinJ = (gMinJ < 0) ? 0 : gMinJ;
	int gMinK {tk - deltaK};
	gMinK = (gMinK < 0) ? 0 : gMinK;

	int gMaxI {ti + deltaI};
	gMaxI = (gMaxI >= nx_) ? nx_ - 1 : gMaxI;
	int gMaxJ {tj + deltaJ};
	gMaxJ = (gMaxJ >= ny_) ? ny_ - 1 : gMaxJ;
	int gMaxK {tk + deltaK};
	gMaxK = (gMaxK >= nz_) ? nz_ - 1 : gMaxK;

	double myRadius {0};
	uint32_t currIndex {0};
	double updatedMax{0};

	for (int ii = gMinI; ii <= gMaxI; ++ii)
		for(int jj = gMinJ; jj <= gMaxJ; ++jj)
			for(int kk = gMinK; kk <= gMaxK; ++kk){
				//check that its nonzero
				currIndex = ArrayIndicesX(nx_,ny_, nz_,
						ii,jj,kk, currIndex);
				if (probArr[currIndex] > 0.) {
					myRadius = calcRadius(ti,tj,tk, ii,jj,kk, dx_, dy_,dz_);
					if (myRadius < radius){
						//just using a quadratic, by definition will vary between 0 an and 1
						probArr[currIndex] *= pow(myRadius/radius, 2);
						if (probArr[currIndex] > updatedMax) updatedMax = probArr[currIndex];
					}
				}

			}


	return updatedMax;

}

template<typename Z>
inline void node_selection<Z>::updateAndSort(std::vector<double>& probArr, std::vector<bool>& okToAddNode) {


	//now update the imgArr_ with probArr
	uint32_t currIndex{0};
	std::cout << "updating the nodes imgArr_ ...\n";
	for (uint32_t i = 0; i < imgArr_.size(); ++i) {
		currIndex = ArrayIndicesX(nx_, ny_, nz_,
				imgArr_[i].i_,imgArr_[i].j_,imgArr_[i].k_, currIndex);
		imgArr_[i].val_ = probArr[currIndex];
	}

    std::cout << "sorting nodes by prob...\n";

     // initialize original index locations
    std::cout << "generated indices\n";
     std::vector<uint32_t> idx(imgArr_.size());
     for (uint32_t i = 0; i != idx.size(); ++i) idx[i] = i;

     std::cout << "sorting indices\n";
      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
           [this](uint32_t i1, uint32_t i2) {return imgArr_[i1].val_ > imgArr_[i2].val_;});

std::cout << "sorting imgArr_\n";
      //now apply that index to imgArr using indices
      std::sort(imgArr_.begin(), imgArr_.end(),
    		  [&idx, this](size_t i1, size_t i2) {return idx[i1] < idx[i2];});

std::cout << "updating okToAddNode\n";
      //now sort the okToAddNode using indices
      std::sort(okToAddNode.begin(), okToAddNode.end(),
           [&idx](uint32_t i1, uint32_t i2)->bool {return idx[i1] < idx[i2];});

}

#endif /* NODE_SELECTION_H_ */

