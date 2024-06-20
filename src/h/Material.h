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

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:
	unsigned int nset; //!< Number of set

	double E; //!< Young's modulus

public:
	//! Virtual deconstructor
	virtual ~CMaterial(){};

	//!	Read material data from stream Input
	virtual bool Read(ifstream &Input) = 0;

	//!	Write material data to Stream
	virtual void Write(COutputter &output) = 0;
};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:
	double Area; //!< Sectional area of a bar element

public:
	//!	Read material data from stream Input
	virtual bool Read(ifstream &Input);

	//!	Write material data to Stream
	virtual void Write(COutputter &output);
};


class CH8Material : public CMaterial
{
public:
	double Nu; //!< Poisson's ratio of a cube element
public:
	//!	Read material data from stream Input
	virtual bool Read(ifstream &Input);

	//!	Write material data to Stream
	virtual void Write(COutputter &output);
};

//!	Material class for Q4 element
class CQ4Material : public CMaterial
{
public:

	double nu;	//!< POISSON RATIO OF Q4 ELEMENT

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};



// Material class for beam element
class CBeamMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a beam element

	//! Define the principal axis of inertia of the section by giving the normalized vector of the y-axis
	double y_axis[3];	//!< Normalized vector of the y-axis of the section respect to the global coordinate system
						// MUST BE PERPENDICULAR TO THE AXIS OF THE BEAM ELEMENT
	double Iz;		//!< Moment of inertia about z-axis
	double Iy;		//!< Moment of inertia about y-axis
					//!< Iyz must be zero, the local coordinate must be Inertia Principal Axis System
	double Jp;		//!< Torsional constant
	double v;		//!< Poisson's ratio

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};


class CT3Material : public CMaterial
{
public:

	double nu;	//!< POISSON RATIO OF T3 ELEMENT

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};


class CShellMaterial : public CMaterial
{
public:

	double nu;			//!< Poisson's ratio

	double thickness;	//!< Thickness of a shell element

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};