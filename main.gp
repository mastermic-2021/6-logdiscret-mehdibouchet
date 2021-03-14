p = 682492462409094395392022581537473179285250139967739310024802121913471471;
g = Mod(6, p);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;



\\ Baby step - Giant step
\\ Génération de la liste un peu trop longue.
\\ N'en parlons pas des prochaines étapes..

\\ g_= g;
\\ B = floor( sqrt(p) );
\\ L = List( Mod(1,p) );
bsgs(A, g, n) = {
    baby= List( g^0 );
    B = floor( sqrt(n) );
    for(i=2, B, listput(baby, baby[i-1]*g) );
    print(baby);
}

\\ On part donc sur une méthode rho de Polard
\\ Renvoie n tq g^n = A
\\ appartenance à l'ensemble G1, G2 ou G3
\\ 
G(i) = {
    n= 0;
    x=binary(i);
    for( j= 1,#x, n+= x[j]);
    lift(Mod(n, 3));
}
f(c) = {
    [xk, a, b]= c;
    Gi= G( lift(xk) );
    if( Gi == 0, xk= xk^2; a*= 2; b*= 2,
        Gi == 1, xk= g*xk; a+= 1,
        Gi == 2, xk= A*xk; b+= 1);
    [xk, a, b];
}
polard(g, A) = {
    x1= g; a1= 1; b1= 0;
    [x2, a2, b2]= f( [x1, a1, b1] );

    while(x1 != x2,
        [x1, a1, b1]= f( [x1, a1, b1] );
        [x2, a2, b2]= f(f( [x2, a2, b2] ));
    );
    (a2- a1) / (b1 - b2);
}

\\ La méthode rho de Polard prend trop de temps à s'executer..
\\ polard(g, A);