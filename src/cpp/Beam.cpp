/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Beam.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBeam::CBeam()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;   // each 3D beam element has 3 displacement components and 3 rotation components at each node
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;

	T_ = new double* [12];
}

//	Desconstructor
CBeam::~CBeam()
{
}

//	Read element data from stream Input
bool CBeam::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBeam::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate beam length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dz^2)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];

	double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);	//	Beam length
	double Lxy = sqrt(DX2[0] + DX2[1]);	//	Beam length projected on x'y' plane
	
// Geometry properties of the beam element
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double k = material_->E * material_->Area / L; //coefficient for stiffness matrix
	double G = material_->E / 2 / (1 + material_->v);



// transform matrix from local to global coordinates
	CalculateTransformMatrix();

//	Calculate element stiffness matrix respect to local coordinate


	// element stiffness matrix in LOCAL COORDINATES
	double ke[12][12] = {0};

	ke[0][0] = k;

	ke[1][1] = 12 * material_->E * material_->Iz / (L * L * L);
	ke[1][5] = 6 * material_->E * material_->Iz / (L * L);
	ke[1][7] = -12 * material_->E * material_->Iz / (L * L * L);
	ke[1][11] = 6 * material_->E * material_->Iz / (L * L);

	ke[2][2] = 12 * material_->E * material_->Iy / (L * L * L);
	ke[2][4] = -6 * material_->E * material_->Iy / (L * L);
	ke[2][8] = -12 * material_->E * material_->Iy / (L * L * L);	
	ke[2][10] = -6 * material_->E * material_->Iy / (L * L);

	ke[3][3] = G * material_->Jp / L;	
	ke[3][9] = -G * material_->Jp / L;

	ke[4][2] = -6 * material_->E * material_->Iy / (L * L);
	ke[4][4] = 4 * material_->E * material_->Iy / L;
	ke[4][8] = 6 * material_->E * material_->Iy / (L * L);
	ke[4][10] = 2 * material_->E * material_->Iy / L;

	ke[5][1] = 6 * material_->E * material_->Iz / (L * L);
	ke[5][5] = 4 * material_->E * material_->Iz / L;
	ke[5][7] = -6 * material_->E * material_->Iz / (L * L);
	ke[5][11] = 2 * material_->E * material_->Iz / L;

	ke[6][0] = -k;
	ke[6][6] = k;

	ke[7][1] = -12 * material_->E * material_->Iz / (L * L * L);
	ke[7][5] = -6 * material_->E * material_->Iz / (L * L);
	ke[7][7] = 12 * material_->E * material_->Iz / (L * L * L);
	ke[7][11] = -6 * material_->E * material_->Iz / (L * L);

	ke[8][2] = -12 * material_->E * material_->Iy / (L * L * L);
	ke[8][4] = 6 * material_->E * material_->Iy / (L * L);
	ke[8][8] = 12 * material_->E * material_->Iy / (L * L * L);
	ke[8][10] = 6 * material_->E * material_->Iy / (L * L);

	ke[9][3] = -G * material_->Jp / L;
	ke[9][9] = G * material_->Jp / L;

	ke[10][2] = -6 * material_->E * material_->Iy / (L * L);
	ke[10][4] = 2 * material_->E * material_->Iy / L;
	ke[10][8] = 6 * material_->E * material_->Iy / (L * L);
	ke[10][10] = 4 * material_->E * material_->Iy / L;

	ke[11][1] = 6 * material_->E * material_->Iz / (L * L);
	ke[11][5] = 2 * material_->E * material_->Iz / L;
	ke[11][7] = -6 * material_->E * material_->Iz / (L * L);
	ke[11][11] = 4 * material_->E * material_->Iz / L;

	// assign the lower triangular part of the element stiffness matrix
	for(unsigned int i = 0; i<12; i++)
	{
		for (unsigned int j = 0; j<12; j++)
		{
			ke[j][i] = ke[i][j];
		}
	}	




    // Global stiffness matrix K_global
    double Ke_global[12][12] = {0};

    // Perform the multiplication T^T * Ke * T to get the global coordinate stiffness matrix
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            for (int k = 0; k < 12; ++k) {
                for (int l = 0; l < 12; ++l) {
                    Ke_global[i][j] += T_[k][i] * ke[k][l] * T_[l][j];
                }
            }
        }
    }
	
	// Store Ke_global to Matrix column by column as it's not a sparse matrix
	unsigned int count = 0;
	for (unsigned int j = 0; j < 12 ; j++)
	{
		for (int i = j; i >= 0; i--)
		{
			Matrix[count++] = Ke_global[i][j];
		}
	}

}

//	Calculate element stress 
void CBeam::ElementStress(double* stress, double* Displacement)
{
	/*
	Description:
		Funtion calculating the normal stress and bending moment of a beam element.
	Variables:
		* stress[0] = N_stress_1		Normal stress at the center of inertia of the section (for linear beam, its constant along x axis)
		* stress[1] = Tx		Torsional moment at gauss point 1
		* stress[2] = My_1		Bending moment about the local y-axis at guass point 1
		* stress[3] = My_2		Bending	moment about the local y-axis at guass point 2
		* stress[4] = Mz_1		Bending moment about the local z-axis at gauss point 1
		* stress[5] = Mz_2		Bending moment about the local z-axis at gauss point 2
		* stress[6] = S_y		Shear force along the local y-axis (for linear beam, shear force distribution is cosntant along x axis)
		* stress[7] = S_z		Shear force along the local z-axis
	*/
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element

	if (sizeof(stress) / sizeof(stress[0]) < 8)
	{
		cout << "Error: stress array for beam stress calculation size is not enough ( <10 )" << endl;
		return;
	}

	// Calculate beam length
	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}
	double L = sqrt(L2);	//	Beam length

	// Calculate local x coordinate for gauss point 1 and 2
	double x_gp[2] = {(-3 / sqrt(3) + 1) * L / 2, (3 / sqrt(3) + 1) * L / 2};

	// Gather local displacements
	double de[12] = {0};
	for (unsigned int i = 0; i < 12; i ++)
	{
		if (LocationMatrix_[i])
			de[i] = Displacement[LocationMatrix_[i]-1];
	} 
	
	// Transfer de to local coordinates (T * de)
	double de_local[12] = {0};
	for (unsigned int i = 0; i < 12; i++)
	{
		for (unsigned int j = 0; j < 12; j++)
		{
			de_local[i] += T_[i][j] * de[j];
		}
	}

	// Calculate stress elements
	stress[0] =  material_->E * (de_local[6] - de_local[0]) / L; // stress due to axial force is constant for linear beam
	
	stress[1] = (material_->Jp * material_->E / 2 / (1 + material_->v) ) * (de_local[9] - de_local[3]) / L; // Torsional moment
	
	stress[2] = ( material_->Iy * material_->E / L2 / L) * 
					( (-12 * x_gp[0] + 6 * L) * de_local[2] + ( 6 * x_gp[0] * L - 4 * L2) * de_local[4] + 
					  ( 12 * x_gp[0] - 6 * L) * de_local[8] + ( 6 * x_gp[0] * L - 2 * L2) * de_local[10] );
	stress[3] = ( material_->Iy * material_->E / L2 / L) * 
					( (-12 * x_gp[1] + 6 * L) * de_local[2] + ( 6 * x_gp[1] * L - 4 * L2) * de_local[4] + 
					  ( 12 * x_gp[1] - 6 * L) * de_local[8] + ( 6 * x_gp[1] * L - 2 * L2) * de_local[10] );

	stress[4] = (-material_->Iz * material_->E / L2 / L) * 
					( (-12 * x_gp[0] + 6 * L) * de_local[1] + (-6 * x_gp[0] * L + 4 * L2) * de_local[5] + 
					  ( 12 * x_gp[0] - 6 * L) * de_local[7] + (-6 * x_gp[0] * L + 2 * L2) * de_local[11] );
	stress[5] = (-material_->Iz * material_->E / L2 / L) * 
					( (-12 * x_gp[1] + 6 * L) * de_local[1] + (-6 * x_gp[1] * L + 4 * L2) * de_local[5] + 
					  ( 12 * x_gp[1] - 6 * L) * de_local[7] + (-6 * x_gp[1] * L + 2 * L2) * de_local[11] );

	stress[6] = (material_->Iz * material_->E / L2 / L) * 
				(-12 * de_local[1] - 6 * L * de_local[5] + 12 * de_local[7] - 6 * L * de_local[11]);
	stress[7] = (material_->Iy * material_->E / L2 / L) * 
				(-12 * de_local[2] + 6 * L * de_local[4] + 12 * de_local[8] - 6 * L * de_local[10]);

}





// Calculate transform matrix
void CBeam::CalculateTransformMatrix()
{
	// Calculate beam length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dz^2)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];

	double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);	//	Beam length

	
// Geometry properties of the beam element
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);	// Pointer to material of the element

// set a local coordinates xyz by the principal y axis of inertia
	/*
		local coordinate xyz is inertia principal axis system
		z axis is generated by the cross product of the principal y axis and the beam axis, 
		it is parallel to the principle z axis of inertia but may not the one because of the origin. 
		But this does not matter when we only consider the angle between z and global axes.
	*/
	double x[3] = { 0 }; // x axis unit vector
	double y[3] = { 0 }; // y axis unit vector
	double z[3] = { 0 }; // z axis unit vector

	// transform matrix from local to global coordinates
	//(12 * 12 matrix, but only store the left upper corner)
	double T_block[3][3]={0};

	T_block[0][0] = x[0] = DX[0] / L;
	T_block[0][1] = x[1] = DX[1] / L;
	T_block[0][2] = x[2] = DX[2] / L;

	T_block[1][0] = y[0] = material_->y_axis[0];
	T_block[1][1] = y[1] = material_->y_axis[1];
	T_block[1][2] = y[2] = material_->y_axis[2];

	// Cross product to get the z axis
	T_block[2][0] = z[0] = x[1] * y[2] - x[2] * y[1];
	T_block[2][1] = z[1] = x[2] * y[0] - x[0] * y[2];
	T_block[2][2] = z[2] = x[0] * y[1] - x[1] * y[0];
	
	// Assemble the global transformation matrix T from T_block
    for (unsigned int blockIdx = 0; blockIdx < 12 / 3; ++blockIdx) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                T_[blockIdx * 3 + i][blockIdx * 3 + j] = T_block[i][j];
            }
        }
    }
}
