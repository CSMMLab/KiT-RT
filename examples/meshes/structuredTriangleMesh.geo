NCells = 10;
//+
Point(1) = {-1, 1, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Point(3) = {1, -1, 0, 1.0};
//+
Point(4) = {-1, -1, 0, 1.0};
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
Transfinite Curve {1, 3} = NCells Using Progression 1;
//+
Transfinite Curve {4, 2} = NCells Using Progression 1;
