/*
 * generate_nodes_from_prob_vertices.cxx
 *
 *  Created on: Mar 4, 2015
 *      Author: jkmikell
 */

#include "node_selection.h"
#include <vector>
#include <boost/assert.hpp>

int main(int argc, char* argv[])
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

	/*
	  uint32_t index = myNodeObj.imgArr_.size()-1;
	  std::cout << "last element: " << myNodeObj.imgArr_[index].val_ << " "
						          << myNodeObj.imgArr_[index].i_ << " "
								  << myNodeObj.imgArr_[index].j_ << " "
								  << myNodeObj.imgArr_[index].k_ << "\n";
    */
	node_selection<double> myProbObj = myNodeObj.getProbability();

	node_selection<double> myVertices = myProbObj.CreateListOfVertices(5.0);
	myVertices.WriteFile(probFilename);

	return 0;
}



