/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "T3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CT3::CT3()
{
	NEN_ = 3;	// Each element has 3 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 18;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CT3::~CT3()
{
}

//	Read element data from stream Input
bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = dynamic_cast<CT3Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];

	return true;
}

//	Write element data to stream
void CT3::Write(COutputter& output)
{
	output << setw(9) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber 
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CT3::ElementStiffness(double* Matrix)
{
    // 清空矩阵，准备存放刚度矩阵
    clear(Matrix, SizeOfStiffnessMatrix());

    // 获取材料参数
    CT3Material* material = dynamic_cast<CT3Material*>(ElementMaterial_); // 获取单元材料的指针
    double E = material->E;  // 弹性模量
    double nu = material->nu; // 泊松比

    // 获取节点坐标 --> Loc
    double Loc[3][3];
    for (unsigned int j = 0; j < 3; j++)
        for (unsigned int i = 0; i < 3; i++)
            Loc[j][i] = nodes_[j]->XYZ[i];

    // 弹性矩阵 D
    double factor = E / (1 - nu * nu);
    double D[3][3] = {
        {factor, nu * factor, 0},
        {nu * factor, factor, 0},
        {0, 0, (1 - nu) * factor / 2}
    };

    // 提取节点坐标
    double X1 = Loc[0][0], Y1 = Loc[0][1];
    double X2 = Loc[1][0], Y2 = Loc[1][1];
    double X3 = Loc[2][0], Y3 = Loc[2][1];

    // 计算 A, B, C
    double A1 = X2 * Y3 - X3 * Y2;
    double A2 = X3 * Y1 - X1 * Y3;
    double A3 = X1 * Y2 - X2 * Y1;

    double B1 = Y2 - Y3, B2 = Y3 - Y1, B3 = Y1 - Y2;

    double C1 = X3 - X2, C2 = X1 - X3, C3 = X2 - X1;

    double Area = 0.5 * (A1 + A2 + A3);

    // 计算B矩阵
    double B[3][6] =
    {
        {B1 / (2 * Area), 0, B2 / (2 * Area), 0, B3 / (2 * Area), 0},
        {0, C1 / (2 * Area), 0, C2 / (2 * Area), 0, C3 / (2 * Area)},
        {C1 / (2 * Area), B1 / (2 * Area), C2 / (2 * Area), B2 / (2 * Area), C3 / (2 * Area), B3 / (2 * Area)}
    };

    // 计算 B 的转置 Bt
    double Bt[6][3] = {0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
            Bt[j][i] = B[i][j];

    // 计算 Bt * D
    double BtD[6][3] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                BtD[i][j] += Bt[i][k] * D[k][j];

    // 计算 BtD * B
    double BtDB[6][6] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j< 6; j++)
            for (int k = 0; k < 3; k++)
                BtDB[i][j] += BtD[i][k] * B[k][j];

    // 将 BtD * B 乘以常数系数并存入 Ke
    double Ke[6][6] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            Ke[i][j] = Area * BtDB[i][j];

    // 打印 K 矩阵
    std::cout << "K matrix:" << std::endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            std::cout << Ke[i][j] << " ";
        }
        std::cout << std::endl;
    }

    //扩展 K 矩阵，新矩阵的第3列，第6列，第9列，第3行，第6行，第9行都是0，其余部分与原矩阵相同
    double Ke_[18][18] = {0};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            Ke_[i][j] = Ke[i][j];
    
    for (int i = 2; i < 4; i++)
        for (int j = 0; j < 2; j++)
            Ke_[i+4][j] = Ke[i][j];
    
    for (int i = 4; i < 6; i++)
        for (int j = 0; j < 2; j++)
            Ke_[i+8][j] = Ke[i][j];
    
    for (int i = 0; i < 2; i++)
        for (int j = 2; j < 4; j++)
            Ke_[i][j+4] = Ke[i][j];

    for (int i = 2; i < 4; i++)
        for (int j = 2; j < 4; j++)
            Ke_[i+4][j+4] = Ke[i][j];
    
    for (int i = 4; i < 6; i++)
        for (int j = 2; j < 4; j++)
            Ke_[i+8][j+4] = Ke[i][j];
    
    for (int i = 0; i < 2; i++)
        for (int j = 4; j < 6; j++)
            Ke_[i][j+8] = Ke[i][j];
    
    for (int i = 2; i < 4; i++)
        for (int j = 4; j < 6; j++)
            Ke_[i+4][j+8] = Ke[i][j];

    for (int i = 4; i < 6; i++)
        for (int j = 4; j < 6; j++)
            Ke_[i+8][j+8] = Ke[i][j];

    // 打印 K 矩阵
    std::cout << "K matrix:" << std::endl;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            std::cout << Ke_[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // 将 Ke 的值存入输出矩阵 Matrix
    int index = 0;
    for(int j = 0; j < 18; j++)
        for(int i = j; i >= 0; i--)
            Matrix[index++] = Ke_[i][j];
}

// Calculate element stress
void CT3::ElementStress(double* stress, double* Displacement)
{
    // 获取材料属性
    CT3Material* material = dynamic_cast<CT3Material*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    // 获取节点坐标 --> Loc
    double Loc[3][3];
    for (unsigned int j = 0; j < 3; j++)
        for (unsigned int i = 0; i < 3; i++)
            Loc[j][i] = nodes_[j]->XYZ[i];

    // 弹性矩阵 D
    double factor = E / (1 - nu * nu);
    double D[3][3] = {
        {factor, nu * factor, 0},
        {nu * factor, factor, 0},
        {0, 0, (1 - nu) * factor / 2}
    };

    // 提取节点坐标
    double X1 = Loc[0][0], Y1 = Loc[0][1];
    double X2 = Loc[1][0], Y2 = Loc[1][1];
    double X3 = Loc[2][0], Y3 = Loc[2][1];

    // 计算 A, B, C
    double A1 = X2 * Y3 - X3 * Y2;
    double A2 = X3 * Y1 - X1 * Y3;
    double A3 = X1 * Y2 - X2 * Y1;

    double B1 = Y2 - Y3, B2 = Y3 - Y1, B3 = Y1 - Y2;

    double C1 = X3 - X2, C2 = X1 - X3, C3 = X2 - X1;

    double Area = 0.5 * (A1 + A2 + A3);

    // 计算B矩阵
    double B[3][6] =
    {
        {B1 / (2 * Area), 0, B2 / (2 * Area), 0, B3 / (2 * Area), 0},
        {0, C1 / (2 * Area), 0, C2 / (2 * Area), 0, C3 / (2 * Area)},
        {C1 / (2 * Area), B1 / (2 * Area), C2 / (2 * Area), B2 / (2 * Area), C3 / (2 * Area), B3 / (2 * Area)}
    };

    // 计算位移
    double de[6];
    for (int j = 0; j < 3; ++j)
    {
        de[2*j] = 0;
        de[2*j+1] = 0;
        if (LocationMatrix_[6*j+1])
            de[2*j+1] = Displacement[LocationMatrix_[6*j+1] - 1];
        if (LocationMatrix_[6*j])
            de[2*j] = Displacement[LocationMatrix_[6*j] - 1];
    }
    
    std::cout << "d array: ";
    for (unsigned int i = 0; i < 6; i++)
        std::cout << de[i] << " ";
    std::cout << std::endl;

    // 计算DB
    double DB[3][6] = {0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 3; k++)
                DB[i][j] += D[i][k] * B[k][j];
    
    // 计算应力
    for (int i = 0; i < 3; i++)
    {
        stress[i] = 0;
        for (int j = 0; j < 6; j++)
            stress[i] += DB[i][j] * de[j];
    }
}