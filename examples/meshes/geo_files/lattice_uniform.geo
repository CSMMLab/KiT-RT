cl_fine = 0.075;
cl_coarse = 0.1;
length =  0.05;
Point(1) = {-3.5, -3.5, 0, cl_coarse};
Point(2) = {3.5, -3.5, 0, cl_coarse};
Point(3) = {-3.5, 3.5, 0, cl_coarse};
Point(4) = {3.5, 3.5, 0, cl_coarse};
Point(5) = {-1.5, -1.5, 0, cl_fine};
Point(6) = {-2.5, -1.5, 0, cl_fine};
Point(7) = {-2.5, -2.5, 0, cl_fine};
Point(8) = {-1.5, -2.5, 0, cl_fine};
Point(9) = {1.5, -1.5, 0, cl_fine};
Point(10) = {1.5, -2.5, 0, cl_fine};
Point(11) = {2.5, -2.5, 0, cl_fine};
Point(12) = {2.5, -1.5, 0, cl_fine};
Point(13) = {-1.5, 1.5, 0, cl_fine};
Point(14) = {-1.5, 0.5, 0, cl_fine};
Point(15) = {-2.5, 0.5, 0, cl_fine};
Point(16) = {-2.5, -0.5, 0, cl_fine};
Point(17) = {-1.5, -0.5, 0, cl_fine};
Point(18) = {-0.5, -1.5, 0, cl_fine};
Point(19) = {-0.5, -2.5, 0, cl_fine};
Point(20) = {0.5, -2.5, 0, cl_fine};
Point(21) = {0.5, -1.5, 0, cl_fine};
Point(22) = {1.5, -0.5, 0, cl_fine};
Point(23) = {2.5, -0.5, 0, cl_fine};
Point(24) = {2.5, 0.5, 0, cl_fine};
Point(25) = {1.5, 0.5, 0, cl_fine};
Point(26) = {1.5, 1.5, 0, cl_fine};
Point(27) = {0.5, 1.5, 0, cl_fine};
Point(28) = {0.5, 0.5, 0, cl_fine};
Point(29) = {-0.5, 0.5, 0, cl_fine};
Point(30) = {-0.5, 1.5, 0, cl_fine};
Point(31) = {-2.5, 1.5, 0, cl_fine};
Point(32) = {-1.5, 2.5, 0, cl_fine};
Point(33) = {-2.5, 2.5, 0, cl_fine};
Point(34) = {2.5, 1.5, 0, cl_fine};
Point(35) = {2.5, 2.5, 0, cl_fine};
Point(36) = {1.5, 2.5, 0, cl_fine};
Point(37) = {-0.5, -0.5, 0, cl_fine};
Point(38) = {0.5, -0.5, 0, cl_fine};
//+
Line(1) = {3, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 4};
//+
Line(4) = {4, 3};
//+
Line(5) = {29, 37};
//+
Line(6) = {37, 38};
//+
Line(7) = {38, 28};
//+
Line(8) = {28, 29};
//+
Line(9) = {29, 14};
//+
Line(10) = {14, 13};
//+
Line(11) = {13, 30};
//+
Line(12) = {30, 29};
//+
Line(13) = {27, 26};
//+
Line(14) = {26, 25};
//+
Line(15) = {25, 28};
//+
Line(16) = {28, 27};
//+
Line(17) = {26, 36};
//+
Line(18) = {36, 35};
//+
Line(19) = {34, 35};
//+
Line(20) = {34, 26};
//+
Line(21) = {24, 23};
//+
Line(22) = {23, 22};
//+
Line(23) = {22, 25};
//+
Line(24) = {25, 24};
//+
Line(25) = {9, 10};
//+
Line(26) = {10, 11};
//+
Line(27) = {11, 12};
//+
Line(28) = {12, 9};
//+
Line(29) = {9, 21};
//+
Line(30) = {21, 38};
//+
Line(31) = {38, 22};
//+
Line(32) = {22, 9};
//+
Line(33) = {20, 21};
//+
Line(34) = {21, 18};
//+
Line(35) = {19, 18};
//+
Line(36) = {19, 20};
//+
Line(37) = {37, 17};
//+
Line(38) = {5, 17};
//+
Line(39) = {5, 18};
//+
Line(40) = {18, 37};
//+
Line(41) = {14, 15};
//+
Line(42) = {15, 16};
//+
Line(43) = {16, 17};
//+
Line(44) = {17, 14};
//+
Line(45) = {32, 13};
//+
Line(46) = {13, 31};
//+
Line(47) = {31, 33};
//+
Line(48) = {33, 32};
//+
Line(49) = {6, 7};
//+
Line(50) = {7, 8};
//+
Line(51) = {8, 5};
//+
Line(52) = {5, 6};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Curve Loop(3) = {8, -12, -11, -10, 41, 42, 43, -38, 39, -35, 36, 33, -29, -32, -22, -21, -24, -14, -13, -16};
//+
Curve Loop(4) = {45, 46, 47, 48};
//+
Curve Loop(5) = {19, -18, -17, -20};
//+
Curve Loop(6) = {28, 25, 26, 27};
//+
Curve Loop(7) = {51, 52, 49, 50};
//+
Plane Surface(2) = {2, 3, 4, 5, 6, 7};
//+
Plane Surface(3) = {4};
//+
Curve Loop(8) = {12, 9, 10, 11};
//+
Plane Surface(4) = {8};
//+
Curve Loop(9) = {44, -9, 5, 37};
//+
Plane Surface(5) = {9};
//+
Curve Loop(10) = {41, 42, 43, 44};
//+
Plane Surface(6) = {10};
//+
Curve Loop(11) = {37, -38, 39, 40};
//+
Plane Surface(7) = {11};
//+
Plane Surface(8) = {7};
//+
Curve Loop(12) = {35, -34, -33, -36};
//+
Plane Surface(9) = {12};
//+
Curve Loop(13) = {30, -6, -40, -34};
//+
Plane Surface(10) = {13};
//+
Curve Loop(14) = {29, 30, 31, 32};
//+
Plane Surface(11) = {14};
//+
Plane Surface(12) = {6};
//+
Curve Loop(15) = {22, 23, 24, 21};
//+
Plane Surface(13) = {15};
//+
Curve Loop(16) = {15, -7, 31, 23};
//+
Plane Surface(14) = {16};
//+
Curve Loop(17) = {13, 14, 15, 16};
//+
Plane Surface(15) = {17};
//+
Plane Surface(16) = {5};
//+
Physical Curve("void", 53) = {1, 2, 3, 4};
