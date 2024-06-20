/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/
#include "H8.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CH8::CH8()
{
    NEN_ = 8; // Each element has 2 nodes
    nodes_ = new CNode *[NEN_];

    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CH8::~CH8()
{
    delete[] nodes_;
    delete[] LocationMatrix_;
}

bool CH8::Read(ifstream &Input, CMaterial *MaterialSets, CNode *NodeList)
{
    unsigned int MSet; // Material property set number
    unsigned int N1, N2, N3, N4, N5, N6, N7, N8;

    Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CH8Material *>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    nodes_[4] = &NodeList[N5 - 1];
    nodes_[5] = &NodeList[N6 - 1];
    nodes_[6] = &NodeList[N7 - 1];
    nodes_[7] = &NodeList[N8 - 1];

    return true;
}
void CH8::Write(COutputter &output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(9) << nodes_[4]->NodeNumber
           << setw(9) << nodes_[5]->NodeNumber
           << setw(9) << nodes_[6]->NodeNumber
           << setw(9) << nodes_[7]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

unsigned int CH8::SizeOfStiffnessMatrix() { return 300; }
// added to increase processing speed

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by column starting from the diagonal element
void CH8::ElementStiffness(double *Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    // material properties
    CH8Material *material_ = dynamic_cast<CH8Material *>(ElementMaterial_); // Pointer to material of the element
    double E = material_->E;
    double nv = material_->Nu;
    double mu = E / (1 + nv);
    double lbd = nv * E / ((1 + nv) * (1 - 2 * nv));

    double gauss[2];
    gauss[0] = -1 / sqrt(3);
    gauss[1] = 1 / sqrt(3);
    int p = 0;

    for (unsigned int m = 0; m < 2; m++)
    {
        for (unsigned int n = 0; n < 2; n++)
        {
            for (unsigned int o = 0; o < 2; o++)
            {
                double xi = gauss[m];
                double eta = gauss[n];
                double zet = gauss[o];
                double G[24];
                double x1 = 1 - xi;
                double y1 = 1 - eta;
                double z1 = 1 - zet;
                double x2 = 1 + xi;
                double y2 = 1 + eta;
                double z2 = 1 + zet;
                G[0] = -0.125 * y1 * z1;
                G[1] = -0.125 * x1 * z1;
                G[2] = -0.125 * x1 * y1;
                G[3] = -G[0];
                G[4] = -0.125 * x2 * z1;
                G[5] = -0.125 * x2 * y1;
                G[6] = 0.125 * y2 * z1;
                G[7] = -G[4];
                G[8] = -0.125 * x2 * y2;
                G[9] = -G[6];
                G[10] = -G[1];
                G[11] = -0.125 * x1 * y2;
                G[12] = -0.125 * y1 * z2;
                G[13] = -0.125 * x1 * z2;
                G[14] = -G[2];
                G[15] = -G[12];
                G[16] = -0.125 * x2 * z2;
                G[17] = -G[5];
                G[18] = 0.125 * y2 * z2;
                G[19] = -G[16];
                G[20] = -G[8];
                G[21] = -G[18];
                G[22] = -G[13];
                G[23] = -G[11];
                // calculate Jacobian at 0,0,0, store row by row
                double J[9] = {0};
                p = 0;
                for (int i = 0; i < 8; i++)
                {
                    double x = nodes_[i]->XYZ[0];
                    double y = nodes_[i]->XYZ[1];
                    double z = nodes_[i]->XYZ[2];
                    for (int j = 0; j < 3; j++)
                    {
                        J[j * 3] += G[p] * x;
                        J[j * 3 + 1] += G[p] * y;
                        J[j * 3 + 2] += G[p++] * z;
                    }
                }

                // calculate inverse of Jacobian, store row by row
                double invJ[9];
                double detJ = J[0] * (J[4] * J[8] - J[5] * J[7]) -
                              J[1] * (J[3] * J[8] - J[5] * J[6]) +
                              J[2] * (J[3] * J[7] - J[4] * J[6]);
                double invDet = 1.0 / detJ;
                invJ[0] = invDet * (J[4] * J[8] - J[5] * J[7]);
                invJ[1] = invDet * (J[2] * J[7] - J[1] * J[8]);
                invJ[2] = invDet * (J[1] * J[5] - J[2] * J[4]);

                invJ[3] = invDet * (J[5] * J[6] - J[3] * J[8]);
                invJ[4] = invDet * (J[0] * J[8] - J[2] * J[6]);
                invJ[5] = invDet * (J[2] * J[3] - J[0] * J[5]);

                invJ[6] = invDet * (J[3] * J[7] - J[4] * J[6]);
                invJ[7] = invDet * (J[1] * J[6] - J[0] * J[7]);
                invJ[8] = invDet * (J[0] * J[4] - J[1] * J[3]);

                // calculate strain matrix B, store column by column
                double B[144] = {0};
                for (int i = 0; i < 24; i += 3)
                {
                    B[i * 6] = B[i * 6 + 9] = B[i * 6 + 16] = invJ[0] * G[i] + invJ[1] * G[i + 1] + invJ[2] * G[i + 2];
                    B[i * 6 + 3] = B[i * 6 + 7] = B[i * 6 + 17] = invJ[3] * G[i] + invJ[4] * G[i + 1] + invJ[5] * G[i + 2];
                    B[i * 6 + 4] = B[i * 6 + 11] = B[i * 6 + 14] = invJ[6] * G[i] + invJ[7] * G[i + 1] + invJ[8] * G[i + 2];
                }

                // calculate stiffness matrix, store column by column
                p = 0;
                for (int j = 0; j < 24; ++j)
                {
                    for (int i = j; i >= 0; i--)
                    {
                        int ib = i * 6 - 1;
                        int jb = j * 6 - 1;
                        Matrix[p++] += (lbd * (B[ib + 1] * (B[jb + 2] + B[jb + 3]) + B[ib + 2] * (B[jb + 1] + B[jb + 3]) + B[ib + 3] * (B[jb + 1] + B[jb + 2])) + (mu + lbd) * (B[ib + 1] * B[jb + 1] + B[ib + 2] * B[jb + 2] + B[ib + 3] * B[jb + 3]) + mu * (B[ib + 4] * B[jb + 4] + B[ib + 5] * B[jb + 5] + B[ib + 6] * B[jb + 6])) * detJ;
                    }
                }
            }
        }
    }
}

//	Calculate element stress
void CH8::ElementStress(double *stress, double *Displacement)
{
    CH8Material *material_ = dynamic_cast<CH8Material *>(ElementMaterial_); // Pointer to material of the element
    double E = material_->E;
    double nv = material_->Nu;
    double mu = E / (1 + nv);
    double lbd = nv * E / ((1 + nv) * (1 - 2 * nv));

    // get element displacement
    double de[24];
    for (int i = 0; i < 24; i++)
    {
        if (LocationMatrix_[i])
        {
            de[i] = Displacement[LocationMatrix_[i] + 1];
        }
        else
        {
            de[i] = 0;
        }
    }
    // calculate gradient of parent element shape functions at gauss point 0,0,0, store column by column
    double G[24];
    G[0] = G[1] = G[2] =
        G[4] = G[5] =
            G[8] =
                G[9] = G[11] =
                    G[12] = G[13] =
                        G[16] =
                            G[21] = -0.125;

    G[3] =
        G[6] = G[7] =
            G[10] =
                G[14] =
                    G[15] = G[17] =
                        G[18] = G[19] = G[20] =
                            G[22] = G[23] = 0.125;
    // calculate J at gauss point, store row by row
    double J[9] = {0};
    int n = 0;
    for (int i = 0; i < 8; i++)
    {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        double z = nodes_[i]->XYZ[2];
        for (int j = 0; j < 3; j++)
        {
            J[j * 3] += G[n] * x;
            J[j * 3 + 1] += G[n] * y;
            J[j * 3 + 2] += G[n++] * z;
        }
    }
    // calculate invJ at gauss point, store row by row
    double invJ[9];
    double detJ = J[0] * (J[4] * J[8] - J[5] * J[7]) -
                  J[1] * (J[3] * J[8] - J[5] * J[6]) +
                  J[2] * (J[3] * J[7] - J[4] * J[6]);
    double invDet = 1.0 / detJ;
    invJ[0] = invDet * (J[4] * J[8] - J[5] * J[7]);
    invJ[1] = invDet * (J[2] * J[7] - J[1] * J[8]);
    invJ[2] = invDet * (J[1] * J[5] - J[2] * J[4]);

    invJ[3] = invDet * (J[5] * J[6] - J[3] * J[8]);
    invJ[4] = invDet * (J[0] * J[8] - J[2] * J[6]);
    invJ[5] = invDet * (J[2] * J[3] - J[0] * J[5]);

    invJ[6] = invDet * (J[3] * J[7] - J[4] * J[6]);
    invJ[7] = invDet * (J[1] * J[6] - J[0] * J[7]);
    invJ[8] = invDet * (J[0] * J[4] - J[1] * J[3]);

    // calculate grad of N at gauss point, store row by row, i.e. N1x,N2x,N3x,...,N8x,N1y,N2y,...,N8z
    double N[24] = {0};
    n = 0;
    for (int i = 0; i < 24; i += 3)
    {
        N[n] = invJ[0] * G[i] + invJ[1] * G[i + 1] + invJ[2] * G[i + 2];
        N[n + 8] = invJ[3] * G[i] + invJ[4] * G[i + 1] + invJ[5] * G[i + 2];
        N[n + 16] = invJ[6] * G[i] + invJ[7] * G[i + 1] + invJ[8] * G[i + 2];
        n++;
    }

    // calculate strain at gauss point
    double eps[6];
    eps[0] = N[0] * de[0] + N[1] * de[3] + N[2] * de[6] + N[3] * de[9] + N[4] * de[12] + N[5] * de[15] + N[6] * de[18] + N[7] * de[21];
    eps[1] = N[8] * de[1] + N[9] * de[4] + N[10] * de[7] + N[11] * de[10] + N[12] * de[13] + N[13] * de[16] + N[14] * de[19] + N[15] * de[22];
    eps[2] = N[16] * de[2] + N[17] * de[5] + N[18] * de[8] + N[19] * de[11] + N[20] * de[14] + N[21] * de[17] + N[22] * de[20] + N[23] * de[23];
    eps[3] = N[0] * de[1] + N[1] * de[4] + N[2] * de[7] + N[3] * de[10] + N[4] * de[13] + N[5] * de[16] + N[6] * de[19] + N[7] * de[22] +
             N[8] * de[0] + N[9] * de[3] + N[10] * de[6] + N[11] * de[9] + N[12] * de[12] + N[13] * de[15] + N[14] * de[18] + N[15] * de[21];
    eps[4] = N[0] * de[2] + N[1] * de[5] + N[2] * de[8] + N[3] * de[11] + N[4] * de[14] + N[5] * de[17] + N[6] * de[20] + N[7] * de[23] +
             N[16] * de[0] + N[17] * de[3] + N[18] * de[6] + N[19] * de[9] + N[20] * de[12] + N[21] * de[15] + N[22] * de[18] + N[23] * de[21];
    eps[5] = N[8] * de[2] + N[9] * de[5] + N[10] * de[8] + N[11] * de[11] + N[12] * de[14] + N[13] * de[17] + N[14] * de[20] + N[15] * de[23] +
             N[16] * de[1] + N[17] * de[4] + N[18] * de[7] + N[19] * de[10] + N[20] * de[13] + N[21] * de[16] + N[22] * de[19] + N[23] * de[22];

    // calculate stress at 2*2*2 gauss points, store point by point; for each point, store by order s11,s22,s33,s12,s13,s23
    stress[0] = (lbd + mu) * eps[0] + lbd * (eps[1] + eps[2]);
    stress[1] = (lbd + mu) * eps[1] + lbd * (eps[0] + eps[2]);
    stress[2] = (lbd + mu) * eps[2] + lbd * (eps[0] + eps[1]);
    stress[3] = 0.5 * mu * eps[3];
    stress[4] = 0.5 * mu * eps[4];
    stress[5] = 0.5 * mu * eps[5];
}
