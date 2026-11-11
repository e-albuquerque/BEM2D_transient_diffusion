// Retangulo 10x8 [-5, 5] x [-4, 4] - Exemplo 2 da Tese (Secao 3.3a, Tabela 3.1)
lc = 0.55;

Point(1) = {-5.0, -4.0, 0.0, lc};
Point(2) = { 5.0, -4.0, 0.0, lc};
Point(3) = { 5.0,  4.0, 0.0, lc};
Point(4) = {-5.0,  4.0, 0.0, lc};
Point(5) = { 0.0,  0.0, 0.0, lc}; // Ponto de monitoramento central A(0,0)

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point{5} In Surface{1};

Physical Curve("bottom") = {1};
Physical Curve("right")  = {2};
Physical Curve("top")    = {3};
Physical Curve("left")   = {4};
Physical Surface("domain") = {1};

Mesh.ElementOrder = 1;
