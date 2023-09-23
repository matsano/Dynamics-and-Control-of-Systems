% Paramètres du système
J = 0.02;
m = 0.6;
sigma = 0.8;
g = 9.81;

% Matrices du système
w = ((m*g^2)/((1+sigma)*J))^(1/4);

A = [0, 1, 0, 0;
    0, 0, -g/(1+sigma), 0;
    0, 0, 0, 1;
    -m*g/(J), 0, 0, 0];

B = [0;
    0;
    0;
    1/J];

R = eye(4);
Q = eye(1);

C = [1, 0, 0, 0;
    0, 0, 1, 0];

% Calcul des valeurs propres du système en boucle ouverte
%K = [0, 0, 0, 0];
%VPs = eig(A - B*K);

% Mappage des valeurs propres
%SYS = ss(A, B, zeros(4), [0; 0; 0; 0]);
%pzmap(SYS)

% Calcul du gain K
%P = [-w, -2*w, -w+1i*w, -w-1i*w];
%K = place(A, B, P);

% Calcul des valeurs propres du système en boucle fermée
%VPs = eig(A - B*K);

% Mappage des valeurs propres
%SYS = ss(A-B*K, B, zeros(4), [0; 0; 0; 0]);
%pzmap(SYS)

% Matrice de commandabilité
%Com = [B, A*B, A*A*B, A*A*A*B];

% Rang de la matrice de commandabilité
%rang = rank(Com);

% Calcul du gain K
%P = [-w, -2*w, -w+1i*w, -w-1i*w];
%K = place(A, B, P);

% Calcul des valeurs propres du système en boucle fermée
%VPs = eig(A - B*K);

% Calcul gain optimal K qui minimise le critère JLQ
K = lqr(A, B, R, Q, zeros);

% Calcul des valeurs propres du système en boucle fermée
%VPs = eig(A - B*K);

% Mappage des valeurs propres
%SYS = ss(A-B*K, B, zeros(4), [0; 0; 0; 0]);
%pzmap(SYS)

% Matrice de observabilité
%O = [C;
%    C*A;
%    C*A*A;
%    C*A*A*A];

% Rang de la matrice de observabilité
%rang = rank(O);

% Calcul du gain L
P = [-w, -3*w, -2*w+1i*w, -2*w-1i*w];
L = place(A', C', P)';

% Calcul des valeurs propres du système en boucle fermée
VPs = eig(A - L*C);
