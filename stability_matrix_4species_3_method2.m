clc;

%% setting model parameters

C = 0.7; %connectance
S = 4; % number of species

% model parameters for mutualistic interactions
r = [-10; 8; 5; -5]; % r_i parameters
s = [0.08; 0.08; 0.08; 0.08]; % s_i parameters
c = [0.08; 0.08; 0.08; 0.08]; % c_i parameters

%% construction of half-normal distribution |N(0,sigma)|
% see file matrix_generator_method2.m

pd = makedist('Normal','mu',0,'sigma',1);

% half-normal distribution |N(0,sigma)|
t = truncate(pd, 0, inf);
z = random(t,S,S);

%% construction of matrix b_ij

b = zeros(S,S);

% According to Allesina and Tang, mutualistic interactions

%  i) For each pair of interactions(bij,bji) i>j,we draw a random value p from a uniform distribution U[0,1];
% ii) If p < or = C, we draw bij and bji independently from a half-normal distribution |N(0,sigma^2)|.
% iii) If p > C,we assign 0 to both bij and bji.
% iv) All diagonal terms,bii,are set to 0.

for i = 2 : S
    for j = 1 : i-1
       p = rand; % Uniform distribution U(0,1)
        if p <= C
            b(i,j) = z(i,j);
            b(j,i) = z(j,i);
        else
            b(i,j) = 0;
            b(j,i) = 0;
        end
     end
end

% all diagonal elements (b_ii) are set to zero

for k = 1 : S
    b(k,k) = 0;
end

% show matrix b_ij

disp('matrix composed of random values drawn from a half-normal distribution |N(0,sigma)|')
disp(z)

disp('the matrix of b_ij elements is')
disp(b)


%% solve the system of differential equations representing the populations

% left-hand side of the system of equation to be solved
% equation #9 in the paper
myfun = @(x) [r(1)-s(1)*x(1)+(1-c(1)*x(1))*(b(1,2)*x(2)+b(1,3)*x(3)+b(1,4)*x(4));
              r(2)-s(2)*x(2)+(1-c(2)*x(2))*(b(2,1)*x(1)+b(2,3)*x(3)+b(2,4)*x(4));
              r(3)-s(3)*x(3)+(1-c(3)*x(3))*(b(3,1)*x(1)+b(3,2)*x(2)+b(3,4)*x(4));
              r(4)-s(4)*x(4)+(1-c(4)*x(4))*(b(4,1)*x(1)+b(4,2)*x(2)+b(4,3)*x(3))];

% make a starting guess at the solution
x0 = [100;100;100;100]; 

% set option to display information after each iteration
options = optimset('Display','iter','MaxIter',400);

% solve the system of equations
[x,fval,exitflag] = fsolve(myfun,x0,options)

%   exitflag - An integer between -4 and 4 indicating why fsolve
%                  stopped. 
%                  Exitflag = 1 indicates success; anything else
%                  indicates that fsolve could not find a solution.

disp('the populations of species at equilibrium are')
disp(x)

%% assistance to solve non-linear equations using Matlab

% https://ocw.mit.edu/courses/mechanical-engineering/2-086-numerical-computation-for-mechanical-engineers-spring-2013...
%.../matlab-exercises/fsolve_MATLAB_Tutorial.m

%% construction of M_ij matrix, Jacobian matrix

M = zeros (S,S);

for i = 1 : S
    for j = 1 : S
        M(i,i) = -x(i)*(s(i)+c(i)*b(i,j)*x(j)); % equation #17 in the paper
        M(i,j) = x(i)*((1-c(i)*x(i))*b(i,j)); % equation #18 in the paper
    end
end

disp('the Jacobian matrix M_ij is')
disp(M)

e = eig(M);
disp('the eigenvalues of the Jacobian matrix M_ij are')
disp(e)

figure
plot(e,'o')
xlabel('Real')
ylabel('Imaginary')
title('Distribution of the eigenvalues of a Jacobian matrix M_{íj}')


%% construction of A_ij matrix

A = zeros (S,S);

for i = 1 : S
    for j = 1 : S
        A(i,i) = -s(i)-c(i)*b(i,j)*x(j); % equation #20 in the paper
        A(i,j) = (1-c(i)*x(i))*b(i,j); % equation #21 in the paper
    end
end

disp('the matrix of A_ij elements is')
disp(A)

%% calculation of parameter d

d = -trace(A)/S; % equation #33 in the paper
disp('the parameter d is')
disp(d)

%% calculation of parameter a

% mean value of the off-diagonal elements of matrix A, 
E1 = (S^2*mean2(A)-trace(A))/((S-1)*S); % E(A_ij) with i different from j
disp('the mean value of the off-diagonal elements of matrix A is')
disp(E1)

% variance of the off-diagonal elements of matrix A, Var(A_ij) with i different from j
A2 = A.*A; % element-wise multiplication in Matlab
E2 = (S^2*mean2(A2)-trace(A2))/((S-1)*S); % E(A_ij^2) with i different from j

E3 = E2 - E1^2; % variance Var(A_ij) = E (A_ij^2)-E(A_ij)^2 with i different from j
disp('the variance of the off-diagonal elements of matrix A is')
disp(E3)

% mean value of the off-diagonal elements of matrix A*A
E4 = 2*((S-1)*S*mean2(A2)-trace(A2))/((S-1)*S); % E(A_ijA_ji) with i different from j
disp('the mean of the off-diagonal elements of matrix A*A is')
disp(E4)

% estimation of parameter a
a = sqrt(S*E3)*(1+E4/E3); % equation #32 in the paper
disp('the parameter a is')
disp(a)

%% test the stability criteria

a_bar = abs(a);
d_bar = abs(d);

if a_bar < d_bar
   disp('the system is stable')
else
    disp('the system is unstable')
end



