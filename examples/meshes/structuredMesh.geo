nCells = 99;
length = 1.0;
//+
Point(1) = {0, length, 0, 1.0};
//+
Point(2) = {length, length, 0, 1.0};
//+
Point(3) = {length, 0, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("void", 5) = {1, 4, 2, 3};

//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Curve {1, 3} = nCells Using Progression 1;
//+
Transfinite Curve {4, 2} = nCells Using Progression 1;
//+
Recombine Surface {1};
