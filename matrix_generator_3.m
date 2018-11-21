clc;

%% setting parameters 

C = 0.7; % connectance
S = 100; % number of species
d = 10; % displacement d

%% construction of distributions 

pd = makedist('Normal','mu',0,'sigma',0.5);
sigma = 0.5;

% half-normal distribution |N(0,sigma)|
t = truncate(pd, 0, inf);
R = random(t,S,S);

%disp('the matrix R is')
%disp(R)
%disp('the mean value of the elements of matrix R is')
%disp(mean2(R))
%disp('the eigenvalues of matrix R are')
%disp(eig(R))

% normal distribution N(0,sigma)
r = random(pd,S,S);

%disp('the matrix r is')
%disp(r)
%disp('the mean value of the elements of matrix r is')
%disp(mean2(r))
%disp('the eigenvalues of matrix r are')
%disp(eig(r))


%% construction of the matrix according to Allesina and Tang

% mutualistic matrix

%  i) For each pair of interactions(bij,bji) i>j,we draw a random value p from a uniform distribution U[0,1];
% ii) If p < or = C, we draw bij and bji independently from a half-normal distribution |N(0,sigma^2)|.
% iii) If p > C,we assign 0 to both bij and bji.
% iv) All diagonal terms,bii,are set to 0.

B = zeros(S,S);

for i = 2:S
    for j = 1:i-1
       p = rand;
        if p <= C
            B(i,j) = R(i,j);
            B(j,i) = R(j,i);
        else
            B(i,j) = 0;
            B(j,i) = 0;
        end
    end
end

% mean value of the the elements of matrix B, 
% los elementos diagonales son cero 
%E1 = mean2(B);
%disp('the average of the elements of matrix B is')
%disp(E1)

% all diagonal elements are set to -d
for k = 1 : S
    B(k,k) = -d;
end

%disp('the mutualistic matrix B is')
%disp(B)
%disp('the mean value of the elements of mutualistic matrix B is')
%disp(mean2(B))

A = eig(B); % autovalores de la matriz B
mean_A = mean(A);
%disp('the eigenvalues of mutualistic matrix B are')
%disp(eig(B))
disp('the average value of the eigenvalues of the mutualistic matrix B is')
disp(mean_A)
% si todo el correcto este número debe ser -d

figure
plot(A,'o'); % representación gráfica de la elipse de autovalores de B
xlabel('Real')
ylabel('Imaginary')
title('Distribution of the eigenvalues of a mutualistic matrix')

% random matrix

%  i) For each pair of interactions(bij,bji) i > j,we draw a random value p from a uniform distribution U[0,1];
% ii) If p < or = C, we draw bij and bji independently from a normal distribution N(0,sigma^2).
% iii) If p > C, we assign 0 to both bij and bji.
% iv) All diagonal terms,bii,are set to 0.

b = zeros(S,S);

for i = 2:S
    for j = 1:i-1
       p = rand;
        if p <= C
            b(i,j) = r(i,j);
            b(j,i) = r(j,i);
        else
            b(i,j) = 0;
            b(j,i) = 0;
        end
    end
end

% mean value of the the elements of matrix B, 
% los elementos diagonales son cero
%E2 = mean2(b);
%disp('the average of the elements of matrix b is')
%disp(E2)
% si todo el correcto, E2 debe ser cero en el límite m -> infinito


% all diagonal elements are set to -d
for k = 1 : S
    b(k,k) = -d;
 end

%disp('the random matrix b is')
%disp(b)
%disp('the mean value of the elements of random matrix b is')
%disp(mean2(b))

a = eig(b); % autovalores de la matriz b
mean_a = mean (a);
%disp('the eigenvalues of the random matrix b are')
%disp(eig(b))
disp('the average value of the eigenvalues of the random matrix b is')
disp(mean_a)
% si todo el correcto este número debe ser -d

figure
plot(a,'o'); % representación gráfica de la elipse de autovalores de b
xlabel('Real')
ylabel('Imaginary')
title('Distribution of the eigenvalues of a random matrix')
% the eigenvalue distribution should be uniform on a circle of radius sigma*sqrt(S*C)

SC = sigma*sqrt(S*C); % normalization factor sigma*sqrt(S*C)

b_rand = b./SC; % normalized random matrix
a_rand = eig(b_rand); % eigenvalues of the normalized random matrix

figure
plot(a_rand,'o'); % representación gráfica de la elipse normalizada de autovalores de b
xlabel('Real')
ylabel('Imaginary')
title('Distribution of the eigenvalues of a normalized random matrix')
% the eigenvalue distribution should be uniform on a unit circle centered at (-d,0) 

%% derivation of the stability criteria

% mutualistic matrix

A_max = max(real(A)); % the right most eigenvalue of the ellipse A

if A_max < 0
   disp('the mutualistic system is numerically stable')
else
    disp('the mutualistic system is numerically unstable')
end

% random matrix

if SC < d
   disp('the random system is stable, according to Allesina criteria')
else
    disp('the random system is unstable, according to Allesina criteria,')
end

a_max = max(real(a)); %the right most eigenvalue of the ellipse a

if a_max < 0
   disp('the random system is numerically stable')
else
    disp('the random system is numerically unstable')
end

