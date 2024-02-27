cl_fine = 0.01;
cl_finer = 0.005;

// Outer points
Point(1) = {0.65, 0.65, 0, cl_fine};
Point(2) = {0., 0.65, 0, cl_fine};
Point(3) = {0.65, 0., 0, cl_fine};
Point(4) = {0., 0., 0, cl_fine};

// Geometry features
// Black
Point(6) = {0.65, 0.6, 0, cl_fine};

Point(13) = {0.65, 0.4, 0, cl_fine};
Point(14) = {0.6, 0.4, 0, cl_fine};
Point(15) = {0.6, 0.0, 0, cl_fine};


// Green (and blue)

Point(19) = {0.2, 0.4, 0, cl_finer};
Point(20) = {0, 0.4, 0, cl_finer};
Point(21) = {0.2, 0, 0, cl_finer};

Point(22) = {0.15, 0.35, 0, cl_finer};
Point(23) = {0.0, 0.35, 0, cl_finer};
Point(24) = {0.15, 0.0, 0, cl_finer};



// Helper points and lines
Point(27) = {0.6, 0.6, 0, cl_fine};



Point(44) = { cl_fine, 0.6 , 0, cl_finer*2};
Point(45) =  { 0, 0.6 , 0, cl_finer*2};
Point(46) = {cl_fine, 0.6 - cl_fine, 0, cl_finer*2};
Point(47) = {0, 0.6 - cl_fine, 0, cl_finer*2};

Point(52) = {0.4 + cl_fine,  cl_fine, 0, cl_finer*2};
Point(53) = {0.4 - cl_fine,  cl_fine, 0, cl_finer*2};
Point(54) = {0.4 + cl_fine,  0, 0, cl_finer*2};
Point(55) = {0.4 - cl_fine,  0, 0, cl_finer*2};

//+
Line(1) = {2, 1};
//+
Line(2) = {1, 6};
//+
Line(3) = {6, 13};
//+
Line(4) = {13, 3};
//+
Line(5) = {3, 15};
//+
Line(6) = {54, 55};
//+
Line(7) = {55, 21};
//+
Line(8) = {21, 24};
//+
Line(9) = {24, 4};
//+
Line(10) = {4, 23};
//+
Line(11) = {23, 20};
//+
Line(12) = {20, 47};
//+
Line(13) = {47, 45};
//+
Line(14) = {45, 2};


Line(15) = {24, 22};
//+
Line(16) = {22, 23};
//+
Line(17) = {20, 19};
//+
Line(18) = {19, 21};
//+
Line(19) = {55, 53};
//+
Line(20) = {53, 52};
//+
Line(21) = {52, 54};
//+
Line(22) = {47, 46};
//+
Line(23) = {46, 44};
//+
Line(24) = {44, 45};
//+
Line(25) = {44, 27};
//+
Line(26) = {27, 6};
//+
Line(27) = {27, 14};
//+
Line(28) = {14, 13};
//+
Line(29) = {54, 15};
//+
Physical Curve("inflow", 60) = {3};
//+
Physical Curve("void", 60) += {1, 2,  4, 5, 6, 7, 8, 9, 10, 11, 12,13,14, 29};
//+//+
Line(30) = {14, 15};
//+
Curve Loop(1) = {10, -16, -15, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {18, 8, 15, 16, 11, 17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -18, -17, 12, 22, 23, 25, 27, 30, -29, -21, -20, -19};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {23, 24, -13, 22};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, 6, 19, 20};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {5, -30, 28, 4};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {28, -3, -26, 27};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {26, -2, -1, -14, -24, 25};
//+
Plane Surface(8) = {8};
