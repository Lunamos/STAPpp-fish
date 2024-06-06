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

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

bool CBeamMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area >> Iz >> Iy >> Jp >> v;	// Young's modulus and section area, moment of inertia, Poisson's ratio

	return true;
}


//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

void CBeamMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << Iz
		   << Iy << setw(16) << Jp << setw(16) << v
		   << endl;
}


