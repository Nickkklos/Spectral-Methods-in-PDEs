% erfc
N  = (1:8).';
t1 = 0.1;
t2 = 0.5;
t3 = 1.0;

J01 = (pi/(4*t1))^0.5*erfc((N)*sqrt(t1));
J05 = (pi/(4*t2))^0.5*erfc((N)*sqrt(t2));
J1  = (pi/(4*t3))^0.5*erfc((N)*sqrt(t3));

T = table(N, J01, J05, J1, ...
    'VariableNames', {'N','J(0.1,N)','J(0.5,N)','J(1,N)'});

disp(T)