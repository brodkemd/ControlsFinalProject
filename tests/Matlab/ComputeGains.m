m = 100.0;
g = 9.81;
Aref = 1413.72;
R = 1.0;

J_xx = 12814.379;
J_yy = 940258.811;
J_zz = 940258.811;
r_E_x = -10.0;

A = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0; 
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0; 
    0, 0, 0, -1, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, -sqrt(2)*g, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, sqrt(2)*g, 0, 0; 
    0, 0, 0, 0, 0, 0, sqrt(2)*g, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, -sqrt(2)/4, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, sqrt(2)/4; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
B = [0, 0, 0; 
    0, 0, 0; 
    0, 0, 0;
    1/m, 0, 0; 
    0, 1/m, 0; 
    0, 0, 1/m; 
    0, 0, 0; 
    0, 0, 0; 
    0, 0, -r_E_x/J_yy; 
    0, r_E_x/J_zz, 0];
C = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0; 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

poles = [-2.121+2.121640638751059j, -2.121-2.121640638751059j, -1.0605+1.0608203193755295j, -1.0605-1.0608203193755295j, -1.7674999999999998+1.7680338656258823j, -1.7674999999999998-1.7680338656258823j, -0.8837499999999999+0.8840169328129411j, -0.8837499999999999-0.8840169328129411j, -1.23725+1.2376237059381177j, -1.23725-1.2376237059381177j]
[K, tol] = place(A, B, poles)
t = C*inv(A-B*K)*B
k_r = (C*inv(A - B * K) * B)

