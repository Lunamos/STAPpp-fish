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