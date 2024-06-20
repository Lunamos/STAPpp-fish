// Shell.h
// created by Wang Qiyao on 2024.6.3

#pragma once

#include "Element.h"

using namespace std;

//! Shell element class
class CShell : public CElement
{
public:

    //!	Constructor
    CShell();

    //!	Desconstructor
    ~CShell();

    //!	Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //!	Write element data to stream
    virtual void Write(COutputter& output);

    //! Calculate transform matrix from local to global coordinate
    void TransformMatrix(double (&T)[24][24]);

    //!	Calculate shape functions and derivatives at a given point (xi, eta)
    void ShapeFunctions(double xi, double eta, double* N, double* dNdxi, double* dNdeta);

    //!	Calculate Jacobian matrix
    void Jacobian(double xi, double eta, double (&J)[2][2]);

    //!	Calculate the B matrix
    void BMatrix_of_plain(double (&J)[2][2],double (&B)[3][8], double xi, double eta);
    void BMatrix_of_plate(double* dNdxi, double* dNdeta, double (&B)[3][12], double xi, double eta);

    //!	Calculate the D matrix
    void DMatrix_of_plain(double E, double mu, double (&D)[3][3]);
    void DMatrix_of_plate(double E, double mu, double t, double (&D)[3][3]);

    //!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //!	Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement);
};