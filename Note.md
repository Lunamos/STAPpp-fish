# Note for reading STAPpp
## Node.h
1. Member of CNode NDF should be inheritaed and modified for BEAM elements. **?For 3D Beam or shell elements, NDF = 5 or 6?** What about B31
2. Remember that after callint DOmain::CalculateEquationNumber(), bcode do not only contain 0s and 1s. Number of global DOF replace all 0s in former bcode matrix and 1s are set to 0.The bcode for entire problem is a list of bcode of each element.
```
# np is the number of node
# dof is the degree of freedom of current node.
NodeList[np].bcode[dof]
```
3. ElementGroup.cpp: Have to add our own Element type in ElementGroup.cpp
In ElementGroup.h,
```
//! Define set of element types
enum ElementTypes
{
    UNDEFINED = 0,
    Bar,    // Bar element
    Q4,     // 4Q element
    T3,     // 3T element
    H8,     // 8H element
    Beam,   // Beam element
    Plate,  // Plate element
    Shell   // Shell elment
};
```
Use these enum as element name.

4. Add material members in Materials.h. Material geometry specs like moment of inertia. 

5. Modify Outputter and add class OutputBeamElements from base class.

6. Add outputter for beam element in Outputter.cpp 
```
void COutputter::OutputBeamElements(unsigned int EleGrp)
{
...
}
void COutputter::OutputElementStress()
{
    ...
    case ElementTypes::Beam: // beam element
    ...
}

```

TODO::

7. Is the orientation and the intertia principle axises should be determined in input file? 

    It could use (0,0,-1) in local coordinate system as default z axis. But in this case, the stifness matrix is including Iyz? 

8. How to calcuate the stress of the Beam? 

    What is included? (Normal stress, 2 bending moment, 2 shearing force, Axial force, shear stress? )


!!!!!! VERY IMPORTANT !!!!!!!
Ip of beam is not Iy^2 + Iz^2, see strucural mechanices page 75
Jp must be input
VECTOR OF PRINCIPAL COORDINATE OF INERTIA OF SECTION MUST BE INPUT


9. Write a input file document.

