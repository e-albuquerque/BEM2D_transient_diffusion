// Disco circular unitario R = 1 centrado em (0, 0) - Exemplo 3 da Tese (Secao 3.4)
lc = 0.14;

Point(1) = { 0.0,  0.0, 0.0, lc}; // Centro r = 0.0
Point(2) = { 1.0,  0.0, 0.0, lc};
Point(3) = { 0.0,  1.0, 0.0, lc};
Point(4) = {-1.0,  0.0, 0.0, lc};
Point(5) = { 0.0, -1.0, 0.0, lc};
Point(6) = { 0.0,  0.8, 0.0, lc}; // Ponto de monitoramento r = 0.8
Point(7) = { 0.0,  0.5, 0.0, lc}; // Ponto de monitoramento r = 0.5

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point{1, 6, 7} In Surface{1};

Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("domain") = {1};

Mesh.ElementOrder = 1;
