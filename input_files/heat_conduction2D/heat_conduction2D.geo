// Gmsh project created on Thu Aug 21 13:16:52 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 100, 0, 1.0};
//+
Point(4) = {0, 100, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Curve {4} = 101 Using Progression 1;
//+
Transfinite Curve {3} = 101 Using Progression 1;
//+
Transfinite Curve {2} = 101 Using Progression 1;
//+
Transfinite Curve {1} = 101 Using Progression 1;
//+
Recombine Surface {1};
