/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! 8H cube element class with reduced integration
class CH8 : public CElement
{
public:
	//!	Constructor
	CH8();

	//!	Desconstructor
	~CH8();

	//!	Read element data from stream Input
	virtual bool Read(ifstream &Input, CMaterial *MaterialSets, CNode *NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter &output);

	void Write(COutputter &output, unsigned int no);

	unsigned int SizeOfStiffnessMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double *Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double *stress, double *Displacement);
};
