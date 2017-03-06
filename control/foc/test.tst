// test the foc
funcprot(0);
exec('foc.sci');

Id_back = 0; 
Iq_back = 0;
theta = 0;
Vdc = 12;
Id = 0; 
Iq = 5;
R = 2;
Ed = 0;
Eq = 0;
ed = 0;
eq = 0;


[Id_back, Iq_back, theta, Ed, ed, Eq, eq] = foc(Id, Iq, Vdc, R, theta, Id_back, Iq_back, Ed, Eq, ed, eq);

