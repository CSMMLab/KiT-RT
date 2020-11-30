// settings Kersting
lc = 0.005;
width = 1.0;
length = 10.0;
start = -2.5;
 //+
Point(1) = {start, width, 0, lc};
//+
Point(2) = {start+length, width, 0, lc};
//+
Point(3) = {start, 0, 0, lc};
//+
Point(4) = {start+length, 0, 0, lc}; 
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("dirichlet") = {4, 2};
//+
Physical Curve("wall_low") = {1};
Physical Curve("wall_up") = {3};
//+
Transfinite Surface {1} = {1, 2, 4, 3};
//+
Transfinite Curve {4, 2} = 1 Using Progression 1;
//+
Transfinite Curve {1, 3} = 1000 Using Progression 1;
//+
Recombine Surface {1};
