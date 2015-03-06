/*
 * generate_nodes_from_prob_vertices.cxx
 *
 *  Created on: Mar 4, 2015
 *      Author: jkmikell
 */


#include "node_selection.h"
#include <iostream>
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

//    node_selection<int16_t> myNodeObj{};
//	myNodeObj.ReadFile(sortedFilename);
//	std::cout << "myNodeObj imgArr size: " << myNodeObj.imgArr_.size() << "\n";
//	node_selection<double> myProbObj = myNodeObj.getProbability();
//	//remove probabilities that are 0
//	uint32_t index{0};
//	std::cout << "there are " << myProbObj.imgArr_.size() << " nodes in myProbObj\n";
//	for (uint32_t i = 0; i < myProbObj.imgArr_.size(); ++i) {
//		if (myProbObj.imgArr_[i].val_ > 0.) index = i;
//	}
//	myProbObj.imgArr_.resize(index+1);
//	std::cout << "there are " << myProbObj.imgArr_.size() << " nodes in myProbObj\n";
//	myProbObj.WriteFile("prob_"+probFilename);

	node_selection<double> myProbObj{};
	myProbObj.ReadFile(sortedFilename);

	node_selection<double> myVertices = myProbObj.CreateListOfVertices(5.0);
	myVertices.WriteFile("inserted_"+probFilename);
	myProbObj.WriteFile("remaining_"+probFilename);

	return 0;
}



