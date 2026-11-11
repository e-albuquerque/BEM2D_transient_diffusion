// Placa quadrada 6x6 [0, 6] x [0, 6] - Exemplo 1 da Tese (Secao 3.2)
lc = 0.45;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {6.0, 0.0, 0.0, lc};
Point(3) = {6.0, 6.0, 0.0, lc};
Point(4) = {0.0, 6.0, 0.0, lc};

// Pontos de monitoramento ao longo da linha media y = 3.0
Point(5) = {1.5, 3.0, 0.0, lc};
Point(6) = {3.0, 3.0, 0.0, lc};
Point(7) = {4.5, 3.0, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point{5, 6, 7} In Surface{1};

Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};
Physical Surface("domain") = {1};

Mesh.ElementOrder = 1;
