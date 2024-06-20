// Shell.cpp
// created by Wang Qiyao on 2024.6.3

#include "Shell.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CShell::CShell()
{
    NEN_ = 4;	// Each element has 4 nodes
    nodes_ = new CNode*[NEN_];
    
    ND_ = 24;   // Each node contributes 6 DOFs
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CShell::~CShell()
{
    delete [] nodes_;
    delete [] ElementMaterial_;
    delete [] LocationMatrix_;
}

//	Read element data from stream Input
bool CShell::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;	// Material property set number
    unsigned int N1, N2, N3, N4;	// Left node number and right node number

    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CShellMaterial*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

    return true;
}

//	Write element data to stream
void CShell::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

//  Calculate transform matrix from local to global coordinate
void CShell::TransformMatrix(double (&T)[24][24])
{
    //  Implement transform matrix for your chosen element type
    //  For example, for a 4-node quadrilateral element:
    double x[4] = {nodes_[0]->XYZ[0], nodes_[1]->XYZ[0], nodes_[2]->XYZ[0], nodes_[3]->XYZ[0]};
    double y[4] = {nodes_[0]->XYZ[1], nodes_[1]->XYZ[1], nodes_[2]->XYZ[1], nodes_[3]->XYZ[1]};
    double z[4] = {nodes_[0]->XYZ[2], nodes_[1]->XYZ[2], nodes_[2]->XYZ[2], nodes_[3]->XYZ[2]};

    double a[3], b[3], c[3];

    a[0] = x[1] - x[0];
    a[1] = y[1] - y[0];
    a[2] = z[1] - z[0];

    b[0] = x[3] - x[0];
    b[1] = y[3] - y[0];
    b[2] = z[3] - z[0];

    double norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= norm;
    a[1] /= norm;
    a[2] /= norm;

    norm = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    b[0] /= norm;
    b[1] /= norm;
    b[2] /= norm;

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];

    // for (int k = 0; k < 4; k++) {
    //     for (int j = 0; j < 4; j++) {
    //         for (int i = 0; i < 3; i++) {
    //             T[i + j * 6][0 + k * 6] = a[i];
    //             T[i + j * 6][1 + k * 6] = b[i];
    //             T[i + j * 6][2 + k * 6] = c[i];

    //             T[i + 3 + j * 6][3 + k * 6] = a[i];
    //             T[i + 3 + j * 6][4 + k * 6] = b[i];
    //             T[i + 3 + j * 6][5 + k * 6] = c[i];
    //         }
    //     }
    // }
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) {
            T[j + 3 * i][0 + 3 * i] += a[j];
            T[j + 3 * i][1 + 3 * i] += b[j];
            T[j + 3 * i][2 + 3 * i] += c[j];
        }
    }
}

//  Calculate shape functions and derivatives at a given point (xi, eta)
void CShell::ShapeFunctions(double xi, double eta, double* N, double* dNdxi, double* dNdeta)
{
    //  Implement shape functions and derivatives for your chosen element type
    //  For example, for a 4-node quadrilateral element:
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);

    dNdxi[0] = -0.25 * (1 - eta);
    dNdxi[1] = 0.25 * (1 - eta);
    dNdxi[2] = 0.25 * (1 + eta);
    dNdxi[3] = -0.25 * (1 + eta);

    dNdeta[0] = -0.25 * (1 - xi);
    dNdeta[1] = -0.25 * (1 + xi);
    dNdeta[2] = 0.25 * (1 + xi);
    dNdeta[3] = 0.25 * (1 - xi);
}

//  Calculate Jacobian matrix at a given point (xi, eta)
void CShell::Jacobian(double xi, double eta, double (&J)[2][2])
{
    //  for a 4-node quadrilateral element:
    double x[4] = {nodes_[0]->XYZ[0], nodes_[1]->XYZ[0], nodes_[2]->XYZ[0], nodes_[3]->XYZ[0]};
    double y[4] = {nodes_[0]->XYZ[1], nodes_[1]->XYZ[1], nodes_[2]->XYZ[1], nodes_[3]->XYZ[1]};
    double z[4] = {nodes_[0]->XYZ[2], nodes_[1]->XYZ[2], nodes_[2]->XYZ[2], nodes_[3]->XYZ[2]};

    double xx = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2) + pow(z[1] - z[0], 2));
    double yy = sqrt(pow(x[3] - x[0], 2) + pow(y[3] - y[0], 2) + pow(z[3] - z[0], 2));

    double a[4] = {0, xx, xx, 0};
    double b[4] = {0, 0, yy, yy};

    // double N[4], dNdxi[4], dNdeta[4];
    // ShapeFunctions(xi, eta, N, dNdxi, dNdeta);

    J[0][0] = xx / 2;
    J[0][1] = 0;
    J[1][0] = 0;
    J[1][1] = yy / 2;

    // for (int i = 0; i < 4; i++) {
    //     J[0][0] += dNdxi[i] * a[i];
    //     J[0][1] += dNdxi[i] * a[i];
    //     J[1][0] += dNdeta[i] * b[i];
    //     J[1][1] += dNdeta[i] * b[i];
    // }
}

//  Calculate B matrix (strain-displacement matrix)
void CShell::BMatrix_of_plain(double (&J)[2][2], double (&B)[3][8], double xi, double eta)
{
    //  4-node quadrilateral element and 3 DOFs per node
    // clear(B, 3 * 8);

    double G[2][4];
    G[0][0] = 0.25 * (eta - 1);
    G[0][1] = 0.25 * (1 - eta);
    G[0][2] = 0.25 * (1 + eta);
    G[0][3] = 0.25 * (-1 - eta);
    G[1][0] = 0.25 * (xi - 1);
    G[1][1] = 0.25 * (-1 - xi);
    G[1][2] = 0.25 * (1 + xi);
    G[1][3] = 0.25 * (1 - xi);

    double J_inv[2][2];
    double det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    J_inv[0][0] = J[1][1] / det_J;
    J_inv[0][1] = -J[0][1] / det_J;
    J_inv[1][0] = -J[1][0] / det_J;
    J_inv[1][1] = J[0][0] / det_J;

    double dNdx[4], dNdy[4];
    // nabla N = J_inv * G
    for (int i = 0; i < 4; i++) {
        dNdx[i] = J_inv[0][0] * G[0][i] + J_inv[0][1] * G[1][i];
        dNdy[i] = J_inv[1][0] * G[0][i] + J_inv[1][1] * G[1][i];
    }

    for (int i = 0; i < 4; i++) {
        B[0][2 * i] = dNdx[i];
        B[1][2 * i + 1] = dNdy[i];
        B[2][2 * i] = dNdy[i];
        B[2][2 * i + 1] = dNdx[i];
    }
}

//  Calculate B matrix (strain-displacement matrix)
void CShell::BMatrix_of_plate(double* dNdxi, double* dNdeta, double (&B)[3][12], double xi, double eta)
{
    //  4-node quadrilateral element and 3 DOFs per node
    // clear(B, 3 * 12);

    double x[4] = {nodes_[0]->XYZ[0], nodes_[1]->XYZ[0], nodes_[2]->XYZ[0], nodes_[3]->XYZ[0]};
    double y[4] = {nodes_[0]->XYZ[1], nodes_[1]->XYZ[1], nodes_[2]->XYZ[1], nodes_[3]->XYZ[1]};
    double z[4] = {nodes_[0]->XYZ[2], nodes_[1]->XYZ[2], nodes_[2]->XYZ[2], nodes_[3]->XYZ[2]};

    double a = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2) + pow(z[1] - z[0], 2)) / 2;
    double b = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2) + pow(z[2] - z[1], 2)) / 2;

    double factor = 1 / (4 * a * b);

    double xiI[4] = {-1, 1, 1, -1};
    double etaI[4] = {-1, -1, 1, 1};
    
    for (int i = 0; i < 4; i++){
        B[0][3 * i] = - factor * 3 * b / a * xiI[i] * xi * (1 + etaI[i] * eta);
        B[0][3 * i + 2] = - factor * b * xiI[i] * (1 + 3 * xiI[i] * xi) * (1 + etaI[i] * eta);
        B[1][3 * i] = - factor * 3 * a / b * etaI[i] * eta * (1 + xiI[i] * xi);
        B[1][3 * i + 1] = factor * a * etaI[i] * (1 + 3 * etaI[i] * eta) * (1 + xiI[i] * xi);
        B[2][3 * i] = factor * xiI[i] * etaI[i] * (4 - 3 * pow(xi, 2) - 3 * pow(eta, 2));
        B[2][3 * i + 1] = factor * b * xiI[i] * (3 * pow(eta, 2) + 2 * etaI[i] * eta - 1);
        B[2][3 * i + 2] = - factor * a * etaI[i] * (3 * pow(xi, 2) + 2 * xiI[i] * xi - 1); 
    }
}

//  Calculate D matrix (constitutive matrix) for plane stress
void CShell::DMatrix_of_plain(double E, double mu, double (&D)[3][3])
{
    //  Calculate D matrix for plane stress
    double factor = E / (1 - mu * mu);
    D[0][0] = factor;
    D[0][1] = factor * mu;
    D[0][2] = 0;
    D[1][0] = factor * mu;
    D[1][1] = factor;
    D[1][2] = 0;
    D[2][0] = 0;
    D[2][1] = 0;
    D[2][2] = factor * (1 - mu) / 2;
}

//  Calculate D matrix (constitutive matrix) for plane stress
void CShell::DMatrix_of_plate(double E, double mu, double t, double (&D)[3][3])
{
    //  Calculate D matrix for plane stress
    double factor = E * t*t*t / (12 * (1 - mu * mu));
    D[0][0] = factor;
    D[0][1] = factor * mu;
    D[0][2] = 0;
    D[1][0] = factor * mu;
    D[1][1] = factor;
    D[1][2] = 0;
    D[2][0] = 0;
    D[2][1] = 0;
    D[2][2] = factor * (1 - mu) / 2;
}

//  Calculate element stiffness matrix 
//  Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CShell::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);
    double E = material_->E;    // Elastic modulus
    double nu = material_->nu;  // Poisson's ratio
    double t = material_->thickness;    // Shell thickness

    // Calculate transform matrix from local to global coordinate
    double T[24][24] = {0};
    TransformMatrix(T);

    /************************************************************************************************************/
    // Calculate stiffness matrix in local coordinate
    // (1) Calculate stiffness matrix for Q4
    double D_plain[3][3];
    DMatrix_of_plain(E, nu, D_plain);
    
    //  Gauss quadrature points and weights for 1x1 integration
    double K_plain[8][8] = {0};
    double K_plate[12][12] = {0};
    
    // double gp[] = {0.0};
    // double gw[] = {2.0};
    double gp[] = {-0.57735027, 0.57735027};
    double gw[] = {1.0, 1.0};

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            //  Gauss point (xi, eta
            double xi = gp[i];
            double eta = gp[j];

            //  Calculate shape functions and derivatives at Gauss point
            double N[4], dNdxi[4], dNdeta[4];
            ShapeFunctions(xi, eta, N, dNdxi, dNdeta);

            //  Calculate Jacobian matrix and its determinant
            double J[2][2];
            Jacobian(xi, eta, J);
            // cout << endl << "J = " << endl
            //         << J[0][0] << " " << J[0][1] << endl
            //         << J[1][0] << " " << J[1][1] << endl;
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            // cout << "detJ = " << detJ << endl;

            // cout << endl << "D of plain " << endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 3; j++) {
            //         cout << D_plain[i][j] << " ";
            //     }
            //     cout << endl;
            // }

            double B_plain[3][8];
            BMatrix_of_plain(J, B_plain, xi, eta);

            // cout << endl << "B of plain " << endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 8; j++) {
            //         cout << B_plain[i][j] << " ";
            //     }
            //     cout << endl;
            // }

            for (int p = 0; p < 8; p++) {
                for (int q = 0; q < 8; q++) {
                    for (int r = 0; r < 3; r++) {
                        for (int s = 0; s < 3; s++) {
                            K_plain[p][q] += B_plain[r][p] * B_plain[s][q] * D_plain[r][s] * detJ * t * gw[i] * gw[j];
                        }
                    }
                }
            }
            
            // (2) Calculate stiffness matrix for Plate
            double D_plate[3][3];
            DMatrix_of_plate(E, nu, t, D_plate);
            //  Gauss quadrature points and weights for 1x1 integration
            
            // double gp = 0.0;
            // double gw = 2.0;

            // double xi = gp;
            // double eta = gp;

            //  Calculate shape functions and derivatives at Gauss point
            // double N[4], dNdxi[4], dNdeta[4];
            ShapeFunctions(xi, eta, N, dNdxi, dNdeta);

            //  Calculate Jacobian matrix and its determinant
            // double J[2][2];
            Jacobian(xi, eta, J);
            detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

            // cout << endl << "D of plate " << endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 3; j++) {
            //         cout << D_plate[i][j] << " ";
            //     }
            //     cout << endl;
            // }

            double B_plate[3][12];
            BMatrix_of_plate(dNdxi, dNdeta, B_plate, xi, eta);

            // cout << endl << "B of plate " << endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 12; j++) {
            //         cout << B_plate[i][j] << " ";
            //     }
            //     cout << endl;
            // }

            for (int p = 0; p < 12; p++) {
                for (int q = 0; q < 12; q++) {
                    for (int r = 0; r < 3; r++) {
                        for (int s = 0; s < 3; s++) {
                            K_plate[p][q] += B_plate[r][p] * B_plate[s][q] * D_plate[r][s] * detJ * gw[i] * gw[j];
                        }
                    }
                }
            }
            // for (int i = 0; i < 12; i++) {
            //     for (int j = 0; j < 12; j++) {
            //         cout << K_plate[i][j] << " ";
            //     }
            //     cout << endl;
            // }

            double test = 0.0;
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    test += B_plate[i][2] * B_plate[j][4] * D_plate[i][j] * detJ * t * gw[i] * gw[j];
                }
            }
            // cout << "test = " << test << endl;
        }
    }

    // (3) Combine K_plain and K_plate to get the stiffness matrix in local coordinate
    double K_local[24][24] = {0.0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            K_local[i * 6 + 0][j * 6 + 0] = K_plain[i * 2 + 0][j * 2 + 0];
            K_local[i * 6 + 0][j * 6 + 1] = K_plain[i * 2 + 0][j * 2 + 1];
            K_local[i * 6 + 1][j * 6 + 0] = K_plain[i * 2 + 1][j * 2 + 0];
            K_local[i * 6 + 1][j * 6 + 1] = K_plain[i * 2 + 1][j * 2 + 1];

            K_local[i * 6 + 2][j * 6 + 2] = K_plate[i * 3 + 0][j * 3 + 0];
            K_local[i * 6 + 2][j * 6 + 3] = K_plate[i * 3 + 0][j * 3 + 1];
            K_local[i * 6 + 2][j * 6 + 4] = K_plate[i * 3 + 0][j * 3 + 2];
            K_local[i * 6 + 3][j * 6 + 2] = K_plate[i * 3 + 1][j * 3 + 0];
            K_local[i * 6 + 3][j * 6 + 3] = K_plate[i * 3 + 1][j * 3 + 1];
            K_local[i * 6 + 3][j * 6 + 4] = K_plate[i * 3 + 1][j * 3 + 2];
            K_local[i * 6 + 4][j * 6 + 2] = K_plate[i * 3 + 2][j * 3 + 0];
            K_local[i * 6 + 4][j * 6 + 3] = K_plate[i * 3 + 2][j * 3 + 1];
            K_local[i * 6 + 4][j * 6 + 4] = K_plate[i * 3 + 2][j * 3 + 2];
        }
    }

    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         cout << K_local[i+12][j] << " ";
    //     }
    //     cout << endl;
    // }

    /************************************************************************************************************/
    // Calculate stiffness matrix in global coordinate
    // K = T * K_local * T^T

    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         cout << K_local[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    double K[24][24] = {0.0};
    for (int p = 0; p < 24; p++) {
        for (int q = 0; q < 24; q++) {
            for (int i = 0; i < 24; i++) {
                for (int j = 0; j < 24; j++) {
                    K[p][q] += T[p][i] * K_local[i][j] * T[q][j];
                }
            }
        }
    }

    // for (int p = 0; p < 24; p++) {
    //     for (int q = 0; q < 24; q++) {
    //         K[p][q] += K_local[p][q];
    //     }
    // }

    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    // cout << endl << "Part of K" << endl;
    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         cout << K[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    for (int i = 0; i < 300; i++) {
        Matrix[i] = 0;
    }

    int itr = 0;
    for (int p = 0; p < 24; p++) {
        for (int q = 0; q <= p; q++) {
            Matrix[itr] = K[p][p - q];
            itr ++;
        }
    }

    // for (int i = 0; i < 12; i++) {
    //     cout << Matrix[i] << " " << endl;
    // }
}

//  Calculate element stress 
void CShell::ElementStress(double* stress, double* Displacement)
{
    CShellMaterial* material_ = dynamic_cast<CShellMaterial*>(ElementMaterial_);
    double E = material_->E;    // Elastic modulus
    double nu = material_->nu;  // Poisson's ratio
    double t = material_->thickness;    // Shell thickness

    // Calculate transform matrix from local to global coordinate
    double T[24][24] = {0};
    TransformMatrix(T);

    double u_global[24] = {0.0};
    for (int i = 0; i < 24; i ++) {
        if (LocationMatrix_[i]) {
            u_global[i] = Displacement[LocationMatrix_[i] - 1];
        }
    }

    // cout << endl << "u_global = " << endl;
    // for (int j = 0; j < 4; j++) {
    //     for (int i = 0; i < 6; i++) {
    //         cout << u_global[i + j * 6] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    // Calculate displacement of gauss points
    double u_local[24] = {0.0}; 
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 24; j++) {
            if (LocationMatrix_[j]) {
                u_local[i] += T[j][i] * Displacement[LocationMatrix_[j] - 1]; // u_local = T^T * u_global
                // u_local[i] += T[j][i] * Displacement[j];
            }
        }
    }
    // cout << endl << "LM = " << endl;
    // c

    // cout << endl << "u_local = " << endl;
    // for (int j = 0; j < 4; j++) {
    //     for (int i = 0; i < 6; i++) {
    //         cout << u_local[i + j * 6] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    // double u_plain[8] = {0};
    // double u_plate[12] = {0};
    double u_plain[8];
    double u_plate[12];

    for (int j = 0; j < 4; j++) {
        u_plain[j * 2 + 0] = u_local[j * 6];
        u_plain[j * 2 + 1] = u_local[j * 6 + 1];
        u_plate[j * 3 + 0] = u_local[j * 6 + 2];
        u_plate[j * 3 + 1] = u_local[j * 6 + 3];
        u_plate[j * 3 + 2] = u_local[j * 6 + 4];
    }

    // cout << endl << "u_plain = " << endl;
    // for (int j = 0; j < 4; j++) {
    //     for (int i = 0; i < 2; i++) {
    //         cout << u_plain[i + j * 2] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    // cout << endl << "u_plate = " << endl;
    // for (int j = 0; j < 4; j++) {
    //     for (int i = 0; i < 3; i++) {
    //         cout << u_plate[i + j * 3] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    //  Calculate Jacobian matrix and its determinant
    double xi = 0, eta = 0;
    double J[2][2] = {0.0};
    //  Calculate shape functions and derivatives at Gauss point
    double N[4], dNdxi[4], dNdeta[4];
    ShapeFunctions(xi, eta, N, dNdxi, dNdeta);

    Jacobian(xi, eta, J);
    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    double B_plain[3][8] = {0.0};
    BMatrix_of_plain(J, B_plain, xi, eta);

    // cout << endl << "B of plain " << endl;
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 8; j++) {
    //         cout << B_plain[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    double B_plate[3][12] = {0.0};
    BMatrix_of_plate(dNdxi, dNdeta, B_plate, xi, eta);

    // Calculate stress in local coordinate
    double strain_plain[3] = {0.0};
    double strain_plate[3] = {0.0};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 8; j++) {
            strain_plain[i] += B_plain[i][j] * u_plain[j];
        }
        for (int j = 0; j < 12; j++) {
            strain_plate[i] += B_plate[i][j] * u_plate[j];
            // cout << j << u_plate[j] << strain_plate[i] << endl;
        }
    }
    // cout << endl << "strain_plain = " << endl;
    // for (int i = 0; i < 3; i++) {
    //     cout << strain_plain[i] << " ";
    // }
    // cout << endl;
    // cout << endl << "strain_plate = " << endl;
    // for (int i = 0; i < 3; i++) {
    //     cout << strain_plate[i] << " ";
    // }
    // cout << endl << endl;

    double stress_local[6] = {0.0};

    double D_plain[3][3];
    DMatrix_of_plain(E, nu, D_plain);
    double D_plate[3][3];
    DMatrix_of_plate(E, nu, t, D_plate);

    for (int i = 0; i < 3; i++) {
        stress_local[0] += D_plain[0][i] * strain_plain[i];
        stress_local[1] += D_plain[1][i] * strain_plain[i];
        // stress_local[5] += D_plain[2][i] * strain_plain[i];
        stress_local[2] += D_plate[0][i] * strain_plate[i];
        stress_local[3] += D_plate[1][i] * strain_plate[i];
        stress_local[4] += D_plate[2][i] * strain_plate[i];
    }

    // Calculate stress in global coordinate
    for (int i = 0; i < 6; i++) {
        stress[i] = 0;
        for (int j = 0; j < 6; j++) {
            stress[i] += T[i][j] * stress_local[j];
        }
    }

    // cout << endl << "T = " << endl;
    // for (int i = 0; i < 12; i++) {
    //     for (int j = 0; j < 12; j++) {
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }
}