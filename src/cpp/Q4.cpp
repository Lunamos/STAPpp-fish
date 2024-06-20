/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 4 nodes, 3 dimension
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ4::~CQ4()
{
}

//	Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	return true;
}

//	Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
			<< setw(9) << nodes_[1]->NodeNumber
			<< setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element4
void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);    // Pointer to material of the element
    double E = material_->E; 
    double nu = material_->nu;     

    double Gauss_Point[2] = {sqrt(3) / 3, -sqrt(3) / 3}; // 高斯点
    double xy[4][3]; // 现在是三维的
    double BtDB[12][12] = {0};
    double K[12][12] = {0};

    // xy矩阵([x,y,z])
    for (unsigned int i = 0; i < 4; i++)
        for (unsigned int j = 0; j < 3; j++)
            xy[i][j] = nodes_[i]->XYZ[j];

    // 计算弹性矩阵D
    double factor = E / (1 - nu * nu);
    double D[3][3] = {0};
    D[0][0] = factor;
    D[0][1] = factor * nu;
    D[1][0] = factor * nu;
    D[1][1] = factor;
    D[2][2] = factor * (1 - nu) / 2;

    // // 打印 D 矩阵
    // std::cout << "D matrix:" << std::endl;
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         std::cout << D[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    for (int s = 0; s < 2; s++) {
        for (int t = 0; t < 2; t++) {
            double eta = Gauss_Point[t];
            double ksi = Gauss_Point[s];
            double G[2][4] = {{eta - 1, 1 - eta, 1 + eta, -1 - eta}, {-1 + ksi, -1 - ksi, 1 + ksi, 1 - ksi}};
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 4; j++) {
                    G[i][j] /= 4.0;
                }
            }
            double J[2][2] = {0};
            double adjJ[2][2] = {0};
            double detJ = 0.0;
            double invJ[2][2] = {0};
            double nablaN[2][4] = {0};
            double B[3][12] = {0}; // 现在是12列

            // 计算雅克比矩阵
            for(int i = 0; i < 2; i++)
                for(int j = 0; j < 2; j++)
                    for(int k = 0; k < 4; k++)
                        J[i][j] += G[i][k] * xy[k][j];

            // // 打印 J 矩阵
            // std::cout << "J matrix:" << std::endl;
            // for (int i = 0; i < 2; i++) {
            //     for (int j = 0; j < 2; j++) {
            //         std::cout << J[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // 计算雅克比矩阵的伴随矩阵和行列式
            adjJ[0][0] = J[1][1];
            adjJ[0][1] = -J[0][1];
            adjJ[1][0] = -J[1][0];
            adjJ[1][1] = J[0][0];
            detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

            // // 打印行列式
            // std::cout << "detJ: " << detJ << std::endl;

            // 计算雅克比矩阵的逆
            for(int i = 0; i < 2; i++)
                for(int j = 0; j < 2; j++)
                    invJ[i][j] = adjJ[i][j] / detJ;

            // // 打印 invJ 矩阵
            // std::cout << "invJ matrix:" << std::endl;
            // for (int i = 0; i < 2; i++) {
            //     for (int j = 0; j < 2; j++) {
            //         std::cout << invJ[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // 计算形函数的梯度
            for(int i = 0; i < 2; i++)
                for(int j = 0; j < 4; j++)
                    for(int k = 0; k < 2; k++)
                        nablaN[i][j] += invJ[i][k] * G[k][j];

            // // 打印 nablaN 矩阵
            // std::cout << "nablaN matrix:" << std::endl;
            // for (int i = 0; i < 2; i++) {
            //     for (int j = 0; j < 4; j++) {
            //         std::cout << nablaN[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // 计算应变矩阵B
            for(int i = 0; i < 4; i++) {
                B[0][3*i] = nablaN[0][i];
                B[1][3*i+1] = nablaN[1][i];
                B[2][3*i] = nablaN[1][i];
                B[2][3*i+1] = nablaN[0][i];
            }

            // // 打印 B 矩阵
            // std::cout << "B matrix:" << std::endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 12; j++) {
            //         std::cout << B[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // 计算BT*D*B并累加到K
            for(int i = 0; i < 12; i++)
                for(int j = 0; j < 12; j++)
                    for(int k = 0; k < 3; k++)
                        for(int l = 0; l < 3; l++)
                            K[i][j] += B[k][i] * D[k][l] * B[l][j] * detJ;
        }
    }

    // // 打印 K 矩阵
    // std::cout << "K matrix:" << std::endl;
    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         std::cout << K[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // 将矩阵的上三角部分存储到一维数组中
    int index = 0;
    for(int j = 0; j < 12; j++)
        for(int i = j; i >= 0; i--)
            Matrix[index++] = K[i][j];
}



//	Calculate element stress 
void CQ4::ElementStress(double* stress, double* Displacement)
{
    CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);    // Pointer to material of the element
    double E = material_->E; 
    double nu = material_->nu;     
    // 计算弹性矩阵D
    double factor = E / (1 - nu * nu);
    double D[3][3] = {0};
    D[0][0] = factor;
    D[0][1] = factor * nu;
    D[1][0] = factor * nu;
    D[1][1] = factor;
    D[2][2] = factor * (1 - nu) / 2;

    double G0[2][4] = {{-0.25, 0.25, 0.25, -0.25}, {-0.25, -0.25, 0.25, 0.25}};
    double xy[4][3]; // 现在是三维的
    double J[2][2] = {0};
    double adjJ[2][2] = {0};
    double detJ = 0.0;
    double invJ[2][2] = {0};
    double nablaN[2][4] = {0};
    double B[3][12] = {0}; // 现在是12列

    // xy矩阵([x,y,z])
    for (unsigned int i = 0; i < 4; i++)
        for (unsigned int j = 0; j < 3; j++)
            xy[i][j] = nodes_[i]->XYZ[j];

    // 计算雅克比矩阵
    for(unsigned int i = 0; i < 2; i++)
        for(unsigned int j = 0; j < 2; j++)
            for(unsigned int k = 0; k < 4; k++)
                J[i][j] += G0[i][k] * xy[k][j];

    // 计算雅克比矩阵的伴随矩阵和行列式
    adjJ[0][0] = J[1][1];
    adjJ[0][1] = -J[0][1];
    adjJ[1][0] = -J[1][0];
    adjJ[1][1] = J[0][0];
    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    // 计算雅克比矩阵的逆
    for(unsigned int i = 0; i < 2; i++)
        for(unsigned int j = 0; j < 2; j++)
            invJ[i][j] = adjJ[i][j] / detJ;

    // 计算形函数的梯度
    for(unsigned int i = 0; i < 2; i++)
        for(unsigned int j = 0; j < 4; j++)
            for(unsigned int k = 0; k < 2; k++)
                nablaN[i][j] += invJ[i][k] * G0[k][j];

    // // 打印 nablaN 矩阵
    // std::cout << "nablaN matrix:" << std::endl;
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         std::cout << nablaN[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // 计算应变矩阵B
    for(unsigned int i = 0; i < 4; i++) {
        B[0][3*i] = nablaN[0][i];
        B[1][3*i+1] = nablaN[1][i];
        B[2][3*i] = nablaN[1][i];
        B[2][3*i+1] = nablaN[0][i];
    }

    double d[12] = {0}; // 位移
    for (unsigned int i = 0; i < 12; i++)
        if (LocationMatrix_[i])
            d[i] = Displacement[LocationMatrix_[i]-1];

    // // 打印d数组
    // std::cout << "d array: ";
    // for (unsigned int i = 0; i < 12; i++)
    //     std::cout << d[i] << " ";
    // std::cout << std::endl;

    // std::cout << "Location Matrix: ";
    // for (unsigned int i = 0; i < 12; i++)
    //     std::cout << LocationMatrix_[i] << " ";
    // std::cout << std::endl;
    // 计算DB = D*B
    double DB[3][12] = {0};
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 12; j++)
            for (unsigned int k = 0; k < 3; k++)
                DB[i][j] += D[i][k] * B[k][j];

    // // 打印DB矩阵
    // std::cout << "DB matrix:" << std::endl;
    // for (unsigned int i = 0; i < 3; i++) {
    //     for (unsigned int j = 0; j < 12; j++)
    //         std::cout << DB[i][j] << " ";
    //     std::cout << std::endl;
    // }

    // 计算stress = DB*d
    for (unsigned int i = 0; i < 3; i++){
        stress[i] = 0.0;
        for (unsigned int j = 0; j < 12; j++){
            stress[i] += DB[i][j] * d[j];
        }
    }

    // // // 打印stress
    // std::cout << "stress: " << stress[0]<< stress[1]<< stress[2] << std::endl;

}

// void CQ4::ElementStress(double* stress, double* Displacement)
// {
// 	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element
// 	double E = material_->E; 
// 	double nu = material_->nu; 	
// 	// 计算弹性矩阵D
// 	double factor = E / (1 - nu * nu);
// 	double D[3][3] = {0};
// 	D[0][0] = factor;
// 	D[0][1] = factor * nu;
// 	D[1][0] = factor * nu;
// 	D[1][1] = factor;
// 	D[2][2] = factor * (1 - nu) / 2;

// 	double G0[2][4] = {{-0.25, 0.25, 0.25, -0.25}, {-0.25, -0.25, 0.25, 0.25}};
// 	double xy[4][2]; 
// 	double J[2][2] = {0};
// 	double adjJ[2][2] = {0};
// 	double detJ;
// 	double invJ[2][2] = {0};
// 	double nablaN[2][4] = {0};
// 	double B[3][8] = {0};
// 	double BtDB[8][8] = {0};
// 	double K[8][8] = {0};

// 	// xy矩阵([x,y])
// 	for (unsigned int i = 0; i < 4; i++)
// 		for (unsigned int j = 0; j < 2; j++)
// 			xy[i][j] = nodes_[i]->XYZ[j] ;

// 	// 计算雅克比矩阵
// 	for(int i = 0; i < 2; i++)
// 		for(int j = 0; j < 2; j++)
// 			for(int k = 0; k < 4; k++)
// 				J[i][j] += G0[i][k] * xy[k][j];

// 	// 计算雅克比矩阵的伴随矩阵和行列式
// 	adjJ[0][0] = J[1][1];
// 	adjJ[0][1] = -J[0][1];
// 	adjJ[1][0] = -J[1][0];
// 	adjJ[1][1] = J[0][0];
// 	detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

// 	// 计算雅克比矩阵的逆
// 	for(int i = 0; i < 2; i++)
// 		for(int j = 0; j < 2; j++)
// 			invJ[i][j] = adjJ[i][j] / detJ;

// 	// 计算形函数的梯度
// 	for(int i = 0; i < 2; i++)
// 		for(int j = 0; j < 4; j++)
// 			for(int k = 0; k < 2; k++)
// 				nablaN[i][j] += invJ[i][k] * G0[k][j];

// 	// 计算应变矩阵B
// 	for(int i = 0; i < 4; i++) {
// 		B[0][2*i] = nablaN[0][i];
// 		B[1][2*i+1] = nablaN[1][i];
// 		B[2][2*i] = nablaN[1][i];
// 		B[2][2*i+1] = nablaN[0][i];
// 	}

// 	double d[8]; // 位移	
// 	for (unsigned int i = 0; i < 8; i++)
// 		d[i] = Displacement[LocationMatrix_[i]-1];

// 	// 计算DB = D*B
// 	double DB[3][8] = {0};
// 	for (int i = 0; i < 3; i++)
// 		for (int j = 0; j < 8; j++)
// 			for (int k = 0; k < 3; k++)
// 				DB[i][j] += D[i][k] * B[k][j];

// 	// 计算stress = DB*d
// 	for (int i = 0; i < 3; i++){
// 		stress[i] = 0.0; // 初始化
// 		for (int j = 0; j < 8; j++)
// 			stress[i] += DB[i][j] * d[j];
// 	}
// }