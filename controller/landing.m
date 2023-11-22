A = [0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,-1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-5.2750165876516446,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,5.2750165876516446,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,5.2750165876516446,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.3535533905932738,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.3535533905932738;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000];
B = [0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0100000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0100000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0100000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000000000000000;
0.0000000000000000,0.0000000000000000,0.0000106353696270;
0.0000000000000000,-0.0000106353696270,0.0000000000000000];
poles = [-0.0212100000000000 + 0.0212164063875106j,
-0.0212100000000000 - 0.0212164063875106j,
-0.0106050000000000 + 0.0106082031937553j,
-0.0106050000000000 - 0.0106082031937553j,
-0.0176750000000000 + 0.0176803386562588j,
-0.0176750000000000 - 0.0176803386562588j,
-0.0088375000000000 + 0.0088401693281294j,
-0.0088375000000000 - 0.0088401693281294j,
-0.0123725000000000 + 0.0123762370593812j,
-0.0123725000000000 - 0.0123762370593812j];
[K, prec] = place(A, B, poles);
save('/home/marekbrodke/Dropbox/School/Flight_Mechanics/ControlsFinalProject/controller/landing.mat', 'K', 'prec');
exit;
