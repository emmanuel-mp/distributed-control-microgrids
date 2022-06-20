
%{
DISTRIBUTED CONTROL OF MICROGRIDS
IMPLEMENTATION OF CONSENSUS ALGORITHM FOR DIFFERENT CONFIGURATIONS OF THE
NETWORK
MAY 2022
BY EMMANUEL AF MOMPREMIER
%}

clear

n = 7; %1000; %number of nodes in the graph considered
t_span = 0:0.0001:15;


% INITIALIZATION OF MATRICES 

A = rand(n,n); %random matrix to host the adjacency matrix 
disp("Random Adjacency Matrix A");
disp(A);

D = zeros(n,n); %zero matrix to host the degree matrix
% disp("Degree Matrix with Zeros D");
% disp(D);

I = eye(n,n); %identity matrix
% disp("Identity Matrix I");
% disp(I);


% CREATE THE ADJACENCY MATRIX 

for i = 1:n
    for j = 1:n


        if A(i,j) > 0.7
            A(i,j) = 1;

        else
            A(i,j) = 0;

        end


    end
end 

disp("ADJACENCY MATRIX A");
disp(A);


% check symmetry of matrix A
A = A + transpose(A);

for i = 1:n
    for j = 1:n


        if A(i,j) == 2
            A(i,j) = 1;

        else
            A(i,j) = A(i,j);

        end

    end

D(i,i) = sum(A(i,:)); % once A is converted into a symmetrical matrix, create the diagonal degree matrix

end 

disp("Symmetrical Adjacency Matrix A");
disp(A);

disp("Degree Matrix D");
disp(D);



% -------- CASE 1 - CONTINUOUS TIME SYSTEM / FIXED TOPOLOGY --------------\
%STANDARD CONSENSUS ALGORITHM

% CREATE LAPLACIAN MATRIX
L = D - A;
disp("Laplacian Matrix L");
disp(L);


%Prove that L is connected (CONSENSUS WILL BE ACHIEVED IF L IS CONNECTED) 
v = eig(L);
disp("Eigenvalue of L");
disp(v);
disp("Second Eigenvalue of L: Lambda_2");
disp(v(2));

% Check whether the network studied is connected
if v(2)  == 0
    disp('This Graph is Not Connected')

else
    disp('YEAH, THE GRAPH IS CONNECTED!!')

end


% SIMULATION OF THE NETWORK
x_init = rand(n,1); %initial value of the imbalance of each node of the network
%t_span = [0,10]; %simulation time

[t,x] = ode23(@(t,x) -L*x, t_span, x_init); %simulation order with the specific solver 

% disp("Imbalances after Simulation");
% disp(x);

%{
ode23 is a three-stage, third-order, Runge-Kutta method. 
ode45 is a six-stage, fifth-order, Runge-Kutta method. 
ode45 does more work per step than ode23, but can take much larger steps. 
For differential equations with smooth solutions, ode45 is often more accurate than ode23
%}

figure(1) %Fig 1 - Continuous time System with Fixed Topology
plot(t,x);
title('CONTINUOUS TIME / STANDARD CONSENSUS');
xlabel('Time');
ylabel('Imbalance of Nodes');




% -------- CASE 2 - CONTINUOUS TIME SYSTEM / FIXED TOPOLOGY --------------
% NORMALIZED CONSENSUS ALGORITHM (DIVIDE BY THE DEGREE OF EACH NODE)
% CONSENSUS ACHIEVED IF NETWORK IS CONNECTED

% CREATE NORMALIZED LAPLACIAN MATRIX
L_norm = D^(-1) * L ;
disp('Normalized Laplacian Matrix L_norm');
disp(L_norm);

%Prove that L_norm is connected
v_norm = eig(L_norm);
disp('Eigenvalue of L_norm');
disp(v_norm);
disp('Second Eigenvalue of L_norm');
disp(v_norm(2));

if v_norm(2)  == 0
    disp('This Graph (NORMALIZED) is Not Connected')

else
    disp('YEAH, THIS GRAPH (NORMALIZED) IS CONNECTED!')

end


% SIMULATION WITH NORMLIZED L_NORM
[t,x] = ode23(@(t,x) -L_norm*x, t_span, x_init);

% disp("Imbalances after Simulation");
% disp(x);

figure(2) %fIG 2 - Normalized Consensus Algo with Fixed Topology with Continuous Time
plot(t,x);
title('CONTINUOUS TIME / NORMALIZED CONSENSUS');
xlabel('Time');
ylabel('Imbalance of Nodes');






% -------- CASE 3 - DISCRETE TIME SYSTEM / FIXED TOPOLOGY --------------
% STANDARD CONSENSUS ALGORITHM
% CONSENSUS ACHIEVED IF NETWORK IS CONNECTED


epsilon = 1/n; %epsilon is sampling time - must be a small value - (slide 59)
x = [];

% TIME VARYING SYSTEM
k = 0; %increment k

%SIMULATION IN DISCRETE TIME
while max(x_init) - min(x_init) >  1e-3 & k < 2000
    
    %generate L
    [L,~] = calculate_L(n);

    x = [x x_init];
    x_init = (I - epsilon*L) * x_init; %x(0+i)

    k = k+1;

     % epsilon = 1 / n (

end

% disp("Imbalances after Simulation");
% disp(x');

figure(3)
plot(x','*');

title('DISCRETE TIME / STANDARD CONSENSUS');
xlabel('Time');
ylabel('Imbalance of Nodes');



% -------- CASE 4 - DISCRETE TIME SYSTEM / FIXED TOPOLOGY --------------
% NORMALIZED CONSENSUS ALGORITHM
% CONSENSUS ACHIEVED IF NETWORK IS CONNECTED

x_init = rand(n,1);

epsilon = 1/n; %epsilon is sampling time - must be a small value - (slide 59)
x_2 = [];

% TIME VARYING SYSTEM
k = 0; %increment k

%SIMULATION IN DISCRETE TIME
while max(x_init) - min(x_init) >  1e-3 & k < 2000
    
    %generate L
    [~,L_norm] = calculate_L(n);

    x_2 = [x_2 x_init];
    x_init = (I - epsilon*L_norm) * x_init; %x(0+i)

    k = k+1;

     % epsilon = 1 / n (

end

% disp("Imbalances after Simulation");
% disp(x_2');

figure(10)
plot(x_2','*');
title('DISCRETE TIME / NORMALIZED CONSENSUS');
xlabel('Time');
ylabel('Imbalance of Nodes');







% -------------- CASE WITH MORE THAN 1 IMBALANCE ---------------------

x_init = rand(n,1);

x=zeros(n,1);
y=zeros(n,1);

%t_span2 = 0:0.01:5;
y_init = rand(n,1);

[t,x] = ode23(@(t,x) -L*x, t_span, x_init);
[t,y] = ode23(@(t,y) -L*y, t, y_init);

figure(4)
plot(x, y, '*');









% -------------- CASE: CREATE A SPECIFIC FORMATION RATHER THAN CONSENSUS ---------------------

x_0 = rand(n,1);
y_0 = rand(n,1);

a = 2.35;
b = 4.23 ;
for k = 1:n

    x_bar(k) = x_0(k) - ((a+b) *cos(k) - b*cos((a/b + 1)*k));   %- cos(2*pi*k/n);
    y_bar(k) = y_0(k)- ((a+b) *sin(k) - b*sin((a/b + 1)*k)); %- sin(2*pi*k/n);

end

[t,X_bar] = ode23(@(t,X_bar) -L*X_bar, t_span, x_bar');
[t,Y_bar] = ode23(@(t,Y_bar) -L*Y_bar, t_span, y_bar');

x=zeros(n,1);
y=zeros(n,1);

for i = 1:n
x(i) = X_bar(end,i) + ((a+b) *cos(i) - b*cos((a/b + 1)*i));
y(i) = Y_bar(end,i) + ((a+b) *sin(i) - b*sin((a/b + 1)*i)); 
end

figure(5)
plot(x,y);


% -------------- CASE: CONTINUOUS DYNAMIC SYSTEM ---------------------


new_n = 3; % analysis of 2 different imbalances at each node (Elec, Gas)

A_char = [0 1 0; -1 0 0; 1 0 0]; %[0 1; -1 0]; %characteristics matrix 
B = [0; 1; 0]; %[0; 1];


x_initial_state = rand(new_n, 1);

x_initial = rand(n*new_n,1); %each node comes with 2 imbalances (5*2)


disp('x_initial');
disp(x_initial);

% transform the synchronization problem into a control/stability problem
% A and B must be controllable to apply the equations of the transformation

C = ctrb(A_char,B);
r = rank(C);

disp('Rank')
disp(r)

%{
In linear algebra, the rank of a matrix A is the dimension of the vector space generated by its columns.
This corresponds to the maximal number of linearly independent columns of A. 
This, in turn, is identical to the dimension of the vector space spanned by its rows.
%}


disp('Length of the initial value matrix')
disp(length(x_initial))


if  r == length(x_initial_state)
disp("A & B ARE CONTROLLABLE!")

else
    disp("Sorry")
end


% Render A  B stable
for i = 1:new_n
D(i,i) = sum(A(i,:));

end


L = D- A;
disp("L");
disp(L);

v = eig(L); %in case of complex eigen values, need to add the conjugate ---- -1-i -1+i (must add this i to prevent 2 equal eigenvalues)
disp('The EigenValues of L');
disp(v);

lambda_2 = v(2);

K = place(A_char, lambda_2*B, [-3, -2, -1]); %the 2 values of K for which 
disp('K');
disp(K);



A_c = kron(eye(n), A_char) - kron(L, B*K);
[t, x] = ode23(@(t,x) A_c*x, t_span, x_initial);


figure(6)
plot(t,x); 

title('CONTINUOUS TIME / STANDARD CONSENSUS / DYNAMIC STATES');
xlabel('Time');
ylabel('Imbalance of Nodes');

%decoupling 
figure(7)
plot(t,x(:, 1:3:end)); 
title('CONTINUOUS TIME / STANDARD CONSENSUS / DYNAMIC STATES');
xlabel('Time');
ylabel('Imbalance of Nodes for State 1');

figure(8)
plot(t,x(:, 2:3:end)); 
title('CONTINUOUS TIME / STANDARD CONSENSUS / DYNAMIC STATES');
xlabel('Time');
ylabel('Imbalance of Nodes for State 2');

figure(9)
plot(t,x(:, 3:3:end)); 
title('CONTINUOUS TIME / STANDARD CONSENSUS / DYNAMIC STATES');
xlabel('Time');
ylabel('Imbalance of Nodes for State 3');



function [L, L_norm] = calculate_L(n)

A = rand(n,n);
D = zeros(n,n);

for i = 1:n
    for j = 1:n


        if A(i,j) > 0.7
            A(i,j) = 1;

        else
            A(i,j) = 0;

        end


    end
end 


% check symmetry of matrix A
A = A + transpose(A);

for i = 1:n
    for j = 1:n


        if A(i,j) == 2
            A(i,j) = 1;

        else
            A(i,j) = A(i,j);

        end

    end

D(i,i) = sum(A(i,:)); 

end 


L = D - A;
L_norm = D^(-1) * L ;


end





