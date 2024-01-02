cl_fine = 0.02;
cl_coarse = 0.075;
cl_boundary = 0.2;
length = 0.075;
Point(1) = {-3.5, -3.5, 0, cl_boundary};
Point(2) = {3.5, -3.5, 0, cl_boundary};
Point(3) = {-3.5, 3.5, 0, cl_boundary};
Point(4) = {3.5, 3.5, 0, cl_boundary};
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
pt_counter = 39;
mid_x = -1;
mid_y = 1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -1;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -1;
mid_y = -1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -1;
mid_y = -2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -2;
mid_y = 2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -2;
mid_y = 1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -2;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -2;
mid_y = -1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = -2;
mid_y = -2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 0;
mid_y = 1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 0;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 0;
mid_y = -1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 0;
mid_y = -2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 1;
mid_y = 1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 1;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 0;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 1;
mid_y = -1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 1;
mid_y = -2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 2;
mid_y = 2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 2;
mid_y = 1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 2;
mid_y = 0;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 2;
mid_y = -1;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
mid_x = 2;
mid_y = -2;
Point(pt_counter) = {mid_x - length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y - length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x - length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Point(pt_counter) = {mid_x + length, mid_y + length, 0, cl_coarse}; // Point to enforce coarser refinement
pt_counter+= 1;
Line(1) = {1, 2};
Line(2) = {3, 1};
Line(3) = {2, 4};
Line(4) = {4, 3};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 5};
Line(18) = {5, 18};
Line(19) = {18, 19};
Line(20) = {19, 20};
Line(21) = {20, 21};
Line(22) = {21, 9};
Line(23) = {9, 22};
Line(24) = {22, 23};
Line(25) = {23, 24};
Line(26) = {24, 25};
Line(27) = {25, 26};
Line(28) = {26, 27};
Line(29) = {27, 28};
Line(30) = {28, 29};
Line(31) = {29, 30};
Line(32) = {30, 13};
Line(33) = {31, 13};
Line(34) = {13, 32};
Line(35) = {32, 33};
Line(36) = {33, 31};
Line(37) = {26, 34};
Line(38) = {34, 35};
Line(39) = {35, 36};
Line(40) = {36, 26};
Line(41) = {14, 29};
Line(42) = {29, 37};
Line(43) = {37, 17};
Line(44) = {17, 14};
Line(45) = {28, 25};
Line(46) = {25, 22};
Line(47) = {22, 38};
Line(48) = {38, 28};
Line(49) = {37, 38};
Line(50) = {38, 21};
Line(51) = {21, 18};
Line(52) = {18, 37};
//+
Physical Curve("void", 61) = {2, 4, 3, 1};
//+
Line(53) = {57, 55};
//+
Line(54) = {56, 55};
//+
Line(55) = {57, 58};
//+
Line(56) = {58, 56};
//+
Line(57) = {61, 62};
//+
Line(58) = {62, 60};
//+
Line(59) = {60, 59};
//+
Line(60) = {59, 61};
//+
Line(61) = {41, 39};
//+
Line(62) = {39, 40};
//+
Line(63) = {40, 42};
//+
Line(64) = {42, 41};
//+
Line(65) = {65, 63};
//+
Line(66) = {63, 64};
//+
Line(67) = {64, 66};
//+
Line(68) = {66, 65};
//+
Line(69) = {45, 43};
//+
Line(70) = {43, 44};
//+
Line(71) = {44, 46};
//+
Line(72) = {46, 45};
//+
Line(73) = {77, 75};
//+
Line(74) = {75, 76};
//+
Line(75) = {76, 78};
//+
Line(76) = {78, 77};
//+
Line(77) = {81, 79};
//+
Line(78) = {79, 80};
//+
Line(79) = {80, 82};
//+
Line(80) = {82, 81};
//+
Line(81) = {85, 83};
//+
Line(82) = {83, 84};
//+
Line(83) = {84, 86};
//+
Line(84) = {86, 85};
//+
Line(85) = {50, 49};
//+
Line(86) = {49, 47};
//+
Line(87) = {47, 48};
//+
Line(88) = {48, 50};
//+
Line(89) = {70, 69};
//+
Line(90) = {69, 67};
//+
Line(91) = {67, 68};
//+
Line(92) = {68, 70};
//+
Line(93) = {93, 91};
//+
Line(94) = {91, 92};
//+
Line(95) = {92, 94};
//+
Line(96) = {94, 93};
//+
Line(97) = {97, 95};
//+
Line(98) = {95, 96};
//+
Line(99) = {96, 98};
//+
Line(100) = {98, 97};
//+
Line(101) = {73, 74};
//+
Line(102) = {74, 72};
//+
Line(103) = {72, 71};
//+
Line(104) = {71, 73};
//+
Line(105) = {53, 54};
//+
Line(106) = {54, 52};
//+
Line(107) = {52, 51};
//+
Line(108) = {51, 53};
//+
Line(109) = {89, 90};
//+
Line(110) = {90, 88};
//+
Line(111) = {88, 87};
//+
Line(112) = {87, 89};
//+
Line(113) = {108, 110};
//+
Line(114) = {110, 109};
//+
Line(115) = {109, 107};
//+
Line(116) = {107, 108};
//+
Line(117) = {103, 104};
//+
Line(118) = {104, 106};
//+
Line(119) = {106, 105};
//+
Line(120) = {105, 103};
//+
Line(121) = {129, 130};
//+
Line(122) = {130, 128};
//+
Line(123) = {128, 127};
//+
Line(124) = {127, 129};
//+
Line(125) = {123, 124};
//+
Line(126) = {124, 126};
//+
Line(127) = {126, 125};
//+
Line(128) = {125, 123};
//+
Line(129) = {122, 121};
//+
Line(130) = {121, 119};
//+
Line(131) = {119, 120};
//+
Line(132) = {120, 122};
//+
Line(133) = {117, 115};
//+
Line(134) = {115, 116};
//+
Line(135) = {116, 118};
//+
Line(136) = {118, 117};
//+
Line(137) = {113, 111};
//+
Line(138) = {111, 112};
//+
Line(139) = {112, 114};
//+
Line(140) = {114, 113};
//+
Curve Loop(1) = {55, 56, 54, -53};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {36, 33, 34, 35};
//+
Plane Surface(2) = {1, 2};
//+
Curve Loop(3) = {57, 58, 59, 60};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {68, 65, 66, 67};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {89, 90, 91, 92};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {101, 102, 103, 104};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {105, 106, 107, 108};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {85, 86, 87, 88};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {72, 69, 70, 71};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {64, 61, 62, 63};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {74, 75, 76, 73};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {80, 77, 78, 79};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {84, 81, 82, 83};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {111, 112, 109, 110};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {115, 116, 113, 114};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {120, 117, 118, 119};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {97, 98, 99, 100};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {96, 93, 94, 95};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {121, 122, 123, 124};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {127, 128, 125, 126};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {129, 130, 131, 132};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {136, 133, 134, 135};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {138, 139, 140, 137};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {32, 13, 41, 31, 32, 13, 41, 31};
//+
Curve Loop(25) = {30, 42, 49, 48};
//+
Plane Surface(25) = {12, 25};
//+
Curve Loop(26) = {31, 32, 13, 41};
//+
Plane Surface(26) = {10, 26};
//+
Curve Loop(27) = {41, 42, 43, 44};
//+
Plane Surface(27) = {9, 27};
//+
Curve Loop(28) = {14, 15, 16, 44};
//+
Plane Surface(28) = {4, 28};
//+
Curve Loop(29) = {43, 17, 18, 52};
//+
Plane Surface(29) = {8, 29};
//+
Curve Loop(30) = {8, 5, 6, 7};
//+
Plane Surface(30) = {6, 30};
//+
Curve Loop(31) = {51, 19, 20, 21};
//+
Plane Surface(31) = {14, 31};
//+
Curve Loop(32) = {51, 52, 49, 50};
//+
Plane Surface(32) = {13, 32};
//+
Curve Loop(33) = {22, 23, 47, 50};
//+
Plane Surface(33) = {16, 33};
//+
Curve Loop(34) = {47, 48, 45, 46};
//+
Plane Surface(34) = {17, 34};
//+
Curve Loop(35) = {45, 27, 28, 29};
//+
Plane Surface(35) = {18, 35};
//+
Curve Loop(36) = {26, 46, 24, 25, 26, 46, 24, 25};
//+

Curve Loop(37) = {39, 40, 37, 38};
//+
Plane Surface(37) = {23, 37};
//+
Curve Loop(38) = {25, 26, 46, 24};
//+
Plane Surface(38) = {21, 38};
//+
Curve Loop(39) = {12, 9, 10, 11};
//+
Plane Surface(39) = {19, 39};
//+
Curve Loop(40) = {4, 2, 1, 3};

//+
Line(141) = {31, 15};
//+
Line(142) = {16, 6};
//+
Line(143) = {8, 19};
//+
Line(144) = {20, 10};
//+
Line(145) = {12, 23};
//+
Line(146) = {24, 34};

//+
Line(147) = {30, 27};
//+
Curve Loop(41) = {141, 15, 142, 6, 7, 143, 20, 144, 10, 11, 145, 25, 146, 38, 39, 40, 28, -147, 32, 34, 35, 36};
//+
Plane Surface(40) = {40, 41};
//+
Curve Loop(42) = {141, -14, -13, -33};
//+
Plane Surface(41) = {3, 42};
//+
Curve Loop(43) = {31, 147, 29, 30};
//+
Plane Surface(42) = {11, 43};
//+
Curve Loop(44) = {37, -146, 26, 27};
//+
Plane Surface(43) = {22, 44};
//+
Curve Loop(45) = {24, -145, 12, 23};
//+
Plane Surface(44) = {20, 45};
//+
Curve Loop(46) = {22, 9, -144, 21};
//+
Plane Surface(45) = {15, 46};
//+
Curve Loop(47) = {18, 19, -143, 8};
//+
Plane Surface(46) = {7, 47};
//+
Curve Loop(48) = {17, 5, -142, 16};
//+
Plane Surface(47) = {5, 48};
//+
Plane Surface(48) = {3, 42};
