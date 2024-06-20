/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

bool CBarMaterial::Read(ifstream &Input)
{
	Input >> nset; // Number of property set

	Input >> E >> Area; // Young's modulus and section area

	return true;
}

void CBarMaterial::Write(COutputter &output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

bool CShellMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> thickness;	// Young's modulus, Poisson ratio and thickness

	return true;
}


void CShellMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << thickness << endl;
}


bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and Poisson ratio

	return true;
}

void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}

bool CT3Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and section area

	return true;
}

void CT3Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}

bool CH8Material::Read(ifstream &Input)
{
	  Input >> nset; // Number of property set

	  Input >> E >> Nu; // Young's modulus and section area
  
  	return true;
}

void CH8Material::Write(COutputter &output)
{
	output << setw(16) << E << setw(16) << Nu << endl;
}

bool CBeamMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area >> Iz >> Iy >> Jp >> v >> y_axis[0] >> y_axis[1] >> y_axis[2];	// Young's modulus and section area, moment of inertia, Poisson's ratio

	return true;
}

void CBeamMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << Iz << setw(16)
		   << Iy << setw(16) << Jp << setw(16) << v << y_axis[0] << y_axis[1] << y_axis[2]
		   << endl;
}
