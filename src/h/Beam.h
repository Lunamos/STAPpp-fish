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

//! Beam element class
class CBeam : public CElement
{
protected:
// Transform matrix between local inertia principal coordinate system and global
	double** T_;

public:

//!	Constructor
	CBeam();

//!	Desconstructor
	~CBeam();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//! Calculate Transform matrix
	void CalculateTransformMatrix();
};
