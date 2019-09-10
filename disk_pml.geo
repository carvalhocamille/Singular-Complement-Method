h =0.1;
H = 5;
angle = Pi/12.;		
C = Cos(angle);	
S = Sin(angle); 
delta = 0.1;
l = 0.3;

Point(1) = {0,0,0,0.1*h};
Point(2) = {1,0,0,h};
Point(3) = {C,S,0,h};
Point(4) = {-1,0,0,h};
Point(5) = {0,-1,0,h};

Point(6) = {l,0,0,h};
Point(7) = {l*C,l*S,0,h};
Point(8) = {-l,0,0,h};
Point(9) = {0,-l,0,h};

// Point(10) = {delta,0,0,h};
// Point(11) = {delta*C,delta*S,0,h};
// Point(12) = {-delta,0,0,h};
// Point(13) = {0,-delta,0,h};

// Line(1) = {1,10};
// Line(10) = {10,6};
Line(100) = {6,2};

// Line(2) = {1,11};
// Line(20) = {11,7};
Line(200) = {7,3};

// disk ext
Circle(3) = {2,1,3};
Circle(4) = {3,1,4};
Circle(5) = {4,1,5};
Circle(6) = {5,1,2};
// disk l
Circle(40) = {6,1,7};
Circle(41) = {7,1,8};
Circle(42) = {8,1,9};
Circle(43) = {9,1,6};
// dosl delta
// Circle(44) = {10,1,11};
// Circle(45) = {11,1,12};
// Circle(46) = {12,1,13};
// Circle(47) = {13,1,10};


// Line Loop(1) = {1,44,-2};
// Line Loop(2) = {10,40,-20,-44};
 Line Loop(3) = {100,3,-200,-40};
// Plane Surface(1) = {1}; // metal
// Plane Surface(2) = {2}; // metal
  Plane Surface(3) = {3}; // metal


// Line Loop(4) = {2,45,46,47,-1};
// Line Loop(5) = {20,41,42,43,-10,-47,-46,-45};
Line Loop(6) = {200,4,5,6,-100,-43,-42,-41};

// Plane Surface(4) = {4}; // vaccum
// Plane Surface(5) = {5}; // vaccum
Plane Surface(6) = {6}; // vaccum

 Transfinite Line {3} = H;
 // Transfinite Line {40} = H;
 // Transfinite Line {44} = H;
 // Transfinite Line {1} = 2*H;
 // Transfinite Line {10} = 2*H;
 Transfinite Line {100} = 2*H;
 // Transfinite Line {20} = 2*H;
 // Transfinite Line {20} = 2*H;
 Transfinite Line {200} = 2*H;
 
 Transfinite Line {4} = 2.9*H;
 Transfinite Line {5} = 2*H;
 Transfinite Line {6} = 2*H;
 Transfinite Line {41} = 2.9*H;
 Transfinite Line {42} = 2*H;
 Transfinite Line {43} = 2*H;
 // Transfinite Line {45} = 2.9*H;
 // Transfinite Line {46} = 2*H;
 // Transfinite Line {47} = 2*H;

Physical Line("Bdy") = {3,4,5,6};
Physical Line("Bdy PML") = {40,41,42,43};
Physical Surface("Metal 0")={3};
// Physical Surface("Metal eta")={2};
// Physical Surface("Metal 1")={1};
Physical Surface("Vaccum 0")={6};
// Physical Surface("Vaccum eta")={5};
// Physical Surface("Vaccum 1")={4};