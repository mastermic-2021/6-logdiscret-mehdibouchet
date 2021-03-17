ALGO_USED= 1;
DEBUG= 0;

if( ALGO_USED == 2, default(parisizemax, 70m));

p = 682492462409094395392022581537473179285250139967739310024802121913471471;
g = Mod(6, p);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;




\\ **********************************************************************************
\\  Infos :                                                                         *
\\                                                                                  *
\\    ALGO_USED: Méthode utilisé pour le calcul de dk dans la Décomposition de P.H  *
\\          =1 - La méthode rho de Pollard logarithmique    (Prend bcp de temps)    *
\\          =2 - La méthode Baby-Step Giant-Step log discret (Foireux)              *
\\          =3 - znlog (Interdit...)                                                *
\\                                                                                  *
\\    DEBUG: Put a 1 pour avoir un debug                                            *
\\ **********************************************************************************



\\ **********************************************************************************
\\  Fonctions de l'exercice :                                                       *
\\                                                                                  *
\\    dlp_bsgs(A,g,n): Baby-Step Giant-Step log discret.                            *
\\      Trouver x tq A=g^x dans un groupe G= <g> d'ordre n                          *
\\                                                                                  *
\\    G(i): détermine l'ensemble d'assignation à un élément i du groupe G           *
\\      Utile pour la mise en pratique de la méthode rho de Pollard                 *
\\                                                                                  *
\\    f(c, g, A): Fonction de marche pseudo-aléatoire de la méthode rho de Pollard  *
\\      Elle détermine x_(k+1) = f(xk)                                              *
\\                                                                                  *
\\    pollard(g, A, p): Fonction de la méthode rho de Pollard logarithmique         *
\\      Trouver x tq A=g^x dans un groupe G= <g> d'ordre p premier                  *
\\                                                                                  *
\\    ph_generic(g, A, n): Décomposition de Pohlig-Hellman                          *
\\      Trouver x tq A=g^x dans un groupe G= <g> d'ordre n                          *
\\      implemente ph(c, A) la décomposition de P.Hellman                           *
\\      pour un groupe d'ordre p^e puissance première                               *
\\                                                                                  *
\\ **********************************************************************************

\\ Dans cet exercice on priviligie l'usage de la Décomposition de Pohlig-Hellman
\\ à celle du rho de Pollard car cette dernière prend trop de temps à s'executer.

dlp_bsgs(A,g,n) = {
  if(A==1,return(0));
  B = sqrtint(n);
  
  baby = vector(B); baby[1] = lift(g^0);
  for(i=2,B, baby[i] = lift(baby[i-1]*g));
  order = vecsort(baby,,1);
  baby = vector(B,i,baby[order[i]]);
  gs = (baby[B]*g)^-1;

  for(a1=0,B-1,
    i = vecsearch(baby,lift(A));
    if(i, return( order[i]+B*a1 ));
    A = A * gs;
    );
  return(-1);
}

G(i) = {
    n= 0;
    x=binary(i);
    for( j= 1,#x, n+= x[j]);
    lift(Mod(n, 3));
}

f(c, g, A) = {
    [xk, a, b]= c;
    Gi= G( lift(xk) );
    if( Gi == 0, xk= xk^2; a*= 2; b*= 2,
        Gi == 1, xk= g*xk; a+= 1,
        Gi == 2, xk= A*xk; b+= 1);
    [xk, a, b];
}
pollard(g, A, p) = {
    x1= 1; a1= 0; b1=0 ;
    [x2, a2, b2]= f( [x1, a1, b1], g, A );
    while(x1 != x2,
        [x1, a1, b1]= f( [x1, a1, b1], g, A );
        [x2, a2, b2]= f( f( [x2, a2, b2], g, A ), g, A);
    );
    if(b1 != b2, lift(Mod( ( (a2- a1) / (b1 - b2) ), p)), 0);
}


ph(c, A)= {
    [g, p, e]= c;
    ak= 0;
    g_= g^(p^(e-1));

    for(k=0, e-1, 
        Ak= (A*g^(-ak))^( p^ ( e-k-1 ) );
        if( ALGO_USED == 1, dk= pollard( g_, Ak, p);,
            ALGO_USED == 2, dk= dlp_bsgs( Ak, g_, p);,
                            dk= znlog(Ak, g_));
        if(dk == -1, print("FAILED: ELEMENT NON TROUVE DANS L"); return(-1););
        ak+= (p^k)*dk;
    );
    if(DEBUG, print(" EGALITY Ak = g^ak ? ", A == g^ak); );
    ak;
}

ph_generic(g, A, n) = {
    factors= factor(n);
    Xi= List();
    for(r=1, matsize(factors)[1], 
        p_i= factors[r, 1];
        e_i= factors[r, 2];
        g_i= g^( n / (p_i ^ e_i) );
        A_i= A^( n / (p_i ^ e_i) );
        x= ph([g_i, p_i, e_i], A_i);

        if(x== -1, return());
        listput( Xi, Mod(x, p_i^e_i ) );
    );
    lift(Mod(lift(chinese(Xi)), n));
}

\\ znlog(A, g) nous renvoie : 
\\ 267758956908173732866643945372404560321057573412753556315673779448856364

A= Mod(A, p);
ph_generic(g, A, p);