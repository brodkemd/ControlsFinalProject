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
poles = [-0.0730410000000000+0.0183058001464017j;
         -0.0730410000000000-0.0183058001464017j;
         -0.0365205000000000+0.0091529000732008j;
         -0.0365205000000000-0.0091529000732008j;
         -0.0608675000000000+0.0152548334553347j;
         -0.0608675000000000-0.0152548334553347j;
         -0.0304337500000000+0.0076274167276674j;
         -0.0304337500000000-0.0076274167276674j;
         -0.0426072500000000+0.0106783834187343j;
         -0.0426072500000000-0.0106783834187343j];
[K, prec] = place(A, B, poles);
save('c:\Users\Delli\Documents\GitHub\ControlsFinalProject\controller\landing.mat', 'K', 'prec');
exit;
