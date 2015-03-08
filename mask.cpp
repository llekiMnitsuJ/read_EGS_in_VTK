/*
 * mask.cpp
 *
 *  Created on: Mar 8, 2015
 *      Author: jkmikell
 */

#include "mask.h"

namespace jkcm {

mask::mask() {
	// TODO Auto-generated constructor stub

}

mask::~mask() {
	// TODO Auto-generated destructor stub
}

void mask::WriteFile(const std::string& filename) const
{
	std::ofstream ofs;
	ofs.open(filename.c_str());
	if(ofs.is_open())
	{
		ofs << name_ << "\n";
		ofs << nxyz_[0] << " " << nxyz_[1] << " " << nxyz_[2] << "\n";

		auto nx = nxyz_[0];
		auto ny = nxyz_[1];
		auto nz = nxyz_[2];
		int32_t currIndex=0;
		//x varies fastest, z slowest
		for (int32_t k = 0; k < nz; ++k)
		{
			for(int32_t j=0; j < ny; ++j)
			{
				for(int32_t i=0; i < nx; ++i)
				{
					currIndex = ArrayIndicesX(nx, ny,nz, i, j, k, currIndex);
					ofs << " " << arr_[currIndex];
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

void mask::ReadFile(const std::string& filename) {
	std::ifstream ifs;
	ifs.open(filename.c_str());
	if (ifs.is_open())
	{
		std::cout << "reading: " << filename << "\n";

		std::string temp;
		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		name_ = temp;
		std::cout << "name: " << name_ << "\n";
		if (!std::getline(ifs, temp)) std::cout << "problem reading " << filename << "\n";
		std::stringstream ss(temp);
		ss >> nxyz_[0] >> nxyz_[1] >> nxyz_[2];
		std::cout << "nx,ny,nz: " << nxyz_[0] << "," << nxyz_[1] << "," << nxyz_[2] << "\n";

		uint32_t size = nxyz_[0]*nxyz_[1]*nxyz_[2];
		arr_.clear();
		arr_.resize(size);
		for(size_t i=0; i < size-1; ++i) ifs >> arr_[i];
	}
	else
	{
		std::cout << "unable to open file " << filename << "\n";
	}

}

} /* namespace jkcm */
