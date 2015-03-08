/*
 * mask.h
 *
 *  Created on: Mar 8, 2015
 *      Author: jkmikell
 */

#ifndef MASK_H_
#define MASK_H_

#include <boost/type.hpp>
#include <string>

namespace jkcm {

class mask {
public:
	mask();
	virtual ~mask();
	std::string name_;
	uint16_t nxyz_[3];
	std::vector<uint16_t> arr_;

	void WriteFile(const std::string& filename) const;
	void ReadFile(const std::string& filename);

	// i,j,k, are zero based.
	// returned 1D index is 0 based
	//stride fastest in x direction (third dimension)
	template<typename T, typename U, typename S>
	S ArrayIndicesX(const U& nx,const U& ny,const U& nz,
			const T& i,const T& j,const T& k,const S& retType){
		S result = i+j*nx+k*nx*ny;
		return result;
	}
};

} /* namespace jkcm */

#endif /* MASK_H_ */
