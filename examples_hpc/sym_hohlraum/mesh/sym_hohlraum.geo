n_recombine = 15;
n_recombine_coarse = 15;

cl_fine = 0.1;

// Outer points
Point(1) = {-0.65, -0.65, 0, cl_fine};
Point(2) = {0.65, -0.65, 0, cl_fine};
Point(3) = {-0.65, 0.65, 0, cl_fine};
Point(4) = {0.65, 0.65, 0, cl_fine};

// Geometry features
// Black
Point(5) = {-0.65, 0.6, 0, cl_fine};
Point(6) = {0.65, 0.6, 0, cl_fine};
Point(7) = {-0.65, -0.6, 0, cl_fine};
Point(8) = {0.65, -0.6, 0, cl_fine};

// Red
Point(9) = {-0.65, 0.4, 0, cl_fine};
Point(10) = {-0.6, 0.4, 0, cl_fine};
Point(11) = {-0.6, -0.4, 0, cl_fine};
Point(12) = {-0.65, -0.4, 0, cl_fine};

Point(13) = {0.65, 0.4, 0, cl_fine};
Point(14) = {0.6, 0.4, 0, cl_fine};
Point(15) = {0.6, -0.4, 0, cl_fine};
Point(16) = {0.65, -0.4, 0, cl_fine};

// Green (and blue)
Point(17) = {-0.2, -0.4, 0, cl_fine};
Point(18) = {-0.2, 0.4, 0, cl_fine};
Point(19) = {0.2, 0.4, 0, cl_fine};
Point(20) = {0.2, -0.4, 0, cl_fine};

Point(21) = {-0.15, -0.35, 0, cl_fine};
Point(22) = {-0.15, 0.35, 0, cl_fine};
Point(23) = {0.15, 0.35, 0, cl_fine};
Point(24) = {0.15, -0.35, 0, cl_fine};




// Lines of basic geometric features
//+
Line(1) = {3, 5};
//+
//+
Line(3) = {6, 4};
//+
Line(5) = {7, 1};
//+
Line(8) = {2, 8};
//+
//+
Line(10) = {12, 11};
//+
//+
Line(13) = {10, 9};
//+
//+
Line(15) = {15, 16};
//+
//+
Line(17) = {13, 14};
//+
//+
Line(19) = {21, 24};
//+
Line(20) = {24, 23};
//+
Line(21) = {23, 22};
//+
Line(22) = {22, 21};
//+
//Line(26) = {19, 18};
//+
Line(28) = {5, 9};
//+
Line(29) = {12, 7};
//+q
Line(30) = {8, 16};
//+
Line(31) = {13, 6};

// Helper points and lines
Point(25) = {-0.6, 0.6, 0, cl_fine};
Point(26) = {-0.6, -0.6, 0, cl_fine};
Point(27) = {0.6, 0.6, 0, cl_fine};
Point(28) = {0.6, -0.6, 0, cl_fine};

Point(29) = {-0.2, 0.6, 0, cl_fine};
Point(30) = {-0.2, -0.6, 0, cl_fine};
Point(31) = {0.2, 0.6, 0, cl_fine};
Point(32) = {0.2, -0.6, 0, cl_fine};

Point(33) = {-0.2, 0.65, 0, cl_fine};
Point(34) = {-0.2, -0.65, 0, cl_fine};
Point(35) = {0.2, 0.65, 0, cl_fine};
Point(36) = {0.2, -0.65, 0, cl_fine};
//+
Point(37) = {-0.6, 0.65, 0, cl_fine};
Point(38) = {-0.6, -0.65, 0, cl_fine};
Point(39) = {0.6, 0.65, 0, cl_fine};
Point(40) = {0.6, -0.65, 0, cl_fine};
//+
Line(32) = {7, 26};
//+
Line(33) = {26, 11};
//+
Line(34) = {11, 10};
//+
Line(35) = {9, 12};
//+
Line(37) = {28, 8};
//+
Line(38) = {28, 15};
//+
Line(40) = {18, 19};
//+
Line(41) = {19, 23};
//+
Line(42) = {19, 20};
//+
Line(43) = {20, 24};
//+
Line(44) = {21, 17};
//+
Line(45) = {17, 20};
//+
Line(46) = {17, 18};
//+
Line(47) = {10, 25};
//+
Line(48) = {25, 5};
//+
Line(50) = {27, 14};
//+
Line(51) = {27, 6};
//+
Line(52) = {13, 16};
//+
Line(53) = {14, 15};
//+
Line(54) = {19, 14};
//+
Line(55) = {20, 15};
//+
Line(56) = {17, 11};
//+
Line(57) = {10, 18};
//+
Line(58) = {25, 29};
//+
Line(59) = {29, 31};
//+
Line(60) = {31, 27};
//+
Line(61) = {22, 18};
//+
Line(62) = {29, 18};
//+
Line(63) = {31, 19};
//+
Line(64) = {17, 30};
//+
Line(65) = {30, 26};
//+
Line(66) = {30, 32};
//+
Line(67) = {32, 28};
//+
Line(68) = {32, 20};
//+
Line(69) = {1, 38};
//+
Line(70) = {38, 26};
//+
Line(71) = {38, 34};
//+
Line(72) = {34, 36};
//+
Line(73) = {36, 40};
//+
Line(74) = {40, 2};
//+
Line(75) = {40, 28};
//+
Line(76) = {36, 32};
//+
Line(77) = {34, 30};
//+
Line(78) = {39, 27};
//+
Line(79) = {39, 4};
//+
Line(80) = {39, 35};
//+
Line(81) = {35, 31};
//+
Line(82) = {35, 33};
//+
Line(83) = {33, 29};
//+
Line(84) = {33, 37};
//+
Line(85) = {37, 25};
//+
Line(86) = {37, 3};
//+
Curve Loop(1) = {48, -1, -86, 85};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {47, 48, 28, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {34, 13, 35, 10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {33, -10, 29, 32};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {70, -32, 5, 69};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {71, 77, 65, -70};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {72, 76, -66, -77};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {73, 75, -67, -76};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {74, 8, -37, -75};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {37, 30, -15, -38};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {15, -52, 17, 53};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {17, -50, 51, -31};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {51, 3, -79, 78};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {60, -78, 80, 81};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {59, -81, 82, 83};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {58, -83, 84, 85};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {57, -62, -58, -47};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {56, 34, 57, -46};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {65, 33, -56, 64};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {66, 68, -45, 64};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {67, 38, -55, -68};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {55, -53, -54, 42};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {54, -50, -60, 63};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {40, -63, -59, 62};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {21, 61, 40, 41};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {46, -61, 22, 44};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {44, 45, 43, -19};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {43, 20, -41, 42};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {19, 20, 21, 22};
//+
Plane Surface(29) = {29};
//+
Physical Curve("void", 87) = {71, 72, 73, 74, 8, 52, 3, 79, 80, 82, 84, 86, 1, 35, 5, 69,29, 30, 31, 28};
//+
// Physical Curve("inflow", 88) = {};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {7};
//+
Transfinite Surface {8};
//+
Transfinite Surface {9};
//+
Transfinite Surface {10};
//+
Transfinite Surface {11};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {14};
//+
Transfinite Surface {15};
//+
Transfinite Surface {16};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {19};
//+
Transfinite Surface {20};
//+
Transfinite Surface {21};
//+
Transfinite Surface {22};
//+
Transfinite Surface {23};
//+
Transfinite Surface {24};
//+
Transfinite Surface {17};
//+
Transfinite Surface {18};
//+
Transfinite Surface {27};
//+
Transfinite Surface {28};
//+
Transfinite Surface {25};
//+
Transfinite Surface {26};
//+
Transfinite Surface {29};
//+
Transfinite Curve {44, 43, 41, 61} = n_recombine Using Progression 1;
// + all vertical
Transfinite Curve {71, 65, 56, 57, 58, 84, 72, 66, 45, 19, 21, 40, 59, 82, 73, 67, 55, 54, 60, 80} = n_recombine_coarse *2 Using Progression 1;
//+  horizontal wide
Transfinite Curve {35, 34, 46, 22, 20, 42, 53, 52} = n_recombine_coarse * 4 Using Progression 1;
//+ horizontal small
Transfinite Curve {28, 47, 62, 63, 50, 31, 30, 38, 68, 64, 33, 29} = n_recombine_coarse Using Progression 1;
//+ inlets
Transfinite Curve {69, 32, 10, 13, 48, 86, 79, 51, 17, 15, 37, 74} = n_recombine / 3 Using Progression 1;
//+ upper and lower bound
Transfinite Curve {5, 70, 77, 76, 75, 8, 3, 78, 81, 83, 85, 1} = n_recombine_coarse / 2 Using Progression 1;

Recombine Surface "*";
// Define meshing options
//Mesh.Algorithm = 6; // Specify meshing algorithm (e.g., Delaunay)
//Mesh.ElementOrder = 1; // Specify element order
// Generate the mesh
//Mesh 2;