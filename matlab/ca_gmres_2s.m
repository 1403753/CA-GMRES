clear;

%{
A = rand(9);
FID = fopen("file.mtx", 'w');
for i = 1:9
    for j = 1:9
        fprintf(FID, "%i %i %e\n", i, j, A(i,j));
    end
end
fclose(FID);
%}
load sparse9x9.mtx;
O = full(spconvert(sparse9x9));

% lambda = [2.514, -0.667+0.092i, -0.667-0.092i];



m = 7;
e_s = [0;0;1];
e1 = [1;0;0];

%edit A(6,4)*(-1)
%%{
A = [
    0.8329    0.3181    0.6473    0.1097    0.7720    0.5254    0.5201;
    0.2564    0.1192    0.5439    0.0636    0.9329    0.5303    0.3477;
    0.6135    0.9398    0.7210    0.4046    0.9727    0.8611    0.1500;
    0.5822    0.6456    0.5225    0.4484    0.1920    0.4849    0.5861;
    0.5407    0.4795    0.9937    0.3658    0.1389    0.3935    0.2621;
    0.8699    0.6393    0.2187    -0.7635    0.6963    0.6714    0.0445;
    0.2648    0.5447    0.1058    0.6279    0.0938    0.7413    0.7549
    ];
%}
%{
A = [
2 0 0 0 0 0 0;
0 5 0 0 1 0 0;
0 0 4 0 0 0 0;
0 0 0 9 0 0 0;
0 0 0 0 10 0 0;
0 0 0 0 0 12 0;
0 0 3 0 0 0 8;
];
%}
%A = rand(m);
x = ones(m,1);
b = A * x;
x0 = zeros(m,1);
%x0 = rand(m,1);
q1_ref = b - A * x0;

[Q_ref, H_ref] = arnoldi(A, q1_ref, m);
eig(H_ref(1:6,1:6))
eigA = eig(H_ref(1:7,1:7));
[leja, outidx] = modified_leja(eigA, m, ones(m,1));
ritz = transpose(leja)
%{
% Ritz values from h1blackletter
ritz(1) = 3.3307;
ritz(2) = -0.5683;
ritz(3) = 0.5878 + 0.2029i;
ritz(4) = 0.5878 - 0.2029i;
ritz(5) = -0.3413 + 0.2466i;
ritz(6) = -0.3413 - 0.2466i;
%}

%{
% Ritz values from h0blackletter
ritz(1) = 3.3296;
ritz(2) = -0.0403;
ritz(3) = 0.5734;
%}

s = 3;
beta = norm(q1_ref);
zeta = [beta; zeros(m-1,1)];
v1 = q1_ref / beta;

B0_ = eye(3);

%{
B0_ = [ritz(1),0,0; B0_];
B0_(2,2) = real(ritz(2));
B0_(3,3) = real(ritz(2));
B0_(2,3) = - imag(ritz(2))^2;

v2 = A*v1 - ritz(1)*v1;
v3 = A*v2 - real(ritz(2))*v2;
v4 = A*v3 - real(ritz(2))*v3 + imag(ritz(2))^2*v2;
%}

%%{
B0_ = [ritz(1),0,0; B0_];
B0_(2,2) = ritz(2);
B0_(3,3) = real(ritz(3));

v2 = A*v1 - ritz(1)*v1;
v3 = A*v2 - ritz(2)*v2;
v4 = A*v3 - real(ritz(3))*v3;
%}

V0_accute = [v2, v3, v4];

V_bad = [v1, A*v1, A*A*v1, A*A*A*v1];
condnmbr_s1 = cond([v1, V0_accute])
bad_condnmbr_s1 = cond(V_bad)

[Q0_, R0_] = qr([v1,V0_accute], 0);

R0_(1,:) = R0_(1,:) * (-1);
Q0_(:,1) = Q0_(:,1) * (-1);

R0 = R0_(1:3,1:3);
Q0 = Q0_(:, 1:3);

Q0_blackletter = Q0_;

h0_blackletter = R0_ * B0_ * inv(R0);
h0blackletter = h0_blackletter(1:3,:);
%verify result
A * Q0;
Q0_* h0_blackletter;
Q_ref * H_ref;

h0 = h0_blackletter(4,3);

v3 = Q0_(:,3);
v4 = Q0_(:,4);

%{ 
%with 6 ritz values + no imaginary overlap
B1_ = eye(3);
B1_ = [ritz(4) 0 0; B1_];
B1_(2,2) = real(ritz(5));
B1_(3,3) = real(ritz(5));
B1_(2,3) = - imag(ritz(5))^2;
B1 = B1_(1:3,1:3);

v5 = A*v4 - ritz(4)*v4;
v6 = A*v5 - real(ritz(5))*v5;
v7 = A*v6 - real(ritz(5))*v6 + imag(ritz(5))^2*v5;
%}

%%{
% with 6 ritz values + imaginary overlap
B1_ = eye(3);
B1_ = [real(ritz(3)) 0 0; B1_];
B1_(2,2) = real(ritz(5));
B1_(3,3) = real(ritz(5));
B1_(2,3) = - imag(ritz(5))^2;
B1 = B1_(1:3,1:3);

v5 = A*(v4) - real(ritz(3))*(v4) + imag(ritz(3))^2*v3;
v6 = A*v5 - real(ritz(5))*v5;
v7 = A*v6 - real(ritz(5))*v6 + imag(ritz(5))^2*v5;
%}

%{ 
%with 3 ritz values
B1_ = eye(3);
B1_ = [ritz(1) 0 0; B1_];
B1_(2,2) = ritz(2);
B1_(3,3) = real(ritz(3));
B1 = B1_(1:3,1:3);

v5 = A*v4 - ritz(1)*(v4);
v6 = A*v5 - ritz(2)*v5;
v7 = A*v6 - real(ritz(3))*v6;
%}

V1_accute = [v5, v6, v7];

v4bad = A*A*A*v1;
V_bad = [v4bad, A*v4bad, A*A*v4bad, A*A*A*v4bad];

condnmbr_s2 = cond([v4, V1_accute])
bad_condnmbr_s2 = cond(V_bad)

R01_blackletter_accute = transpose(Q0_blackletter) * V1_accute;

V1_accute_accute = V1_accute - Q0_blackletter * R01_blackletter_accute;

[Q1_accute, R1_accute] = qr(V1_accute_accute, 0);

R1 = eye(3);
R1(1,:) = [1, R01_blackletter_accute(4,1), R01_blackletter_accute(4,2)];
R1(2:3,2:3) = R1_accute(1:2, 1:2);

R01_blackletter = transpose(Q0) * [v4, V1_accute(:,1:3)];
R01blackletter = transpose(Q0) * [v4, V1_accute(:,1:2)];

h01_blackletter = -h0_blackletter(1:3,1:3) * R01blackletter * inv(R1) + R01_blackletter * B1_ * inv(R1) - imag(ritz(3))^2*e_s*e1'*inv(R1); %e_s should be e_sk 
h01_blackletter
roh1 = R1_accute(3,3);

roh1_tilde_inv = inv(R1(3,3));

z1 = R01_blackletter_accute(4,3);
z1 = [z1;R1_accute(1:2,3)];

Q1_blackletter = [Q0_blackletter, Q1_accute];
H1 = R1 * B1 * inv(R1) + roh1_tilde_inv*1*z1*transpose(e_s) - h0 * e1 * transpose(e_s) * R01blackletter * inv(R1);

h1 = roh1_tilde_inv * roh1 * 1;
h1_blackletter = [h0_blackletter(1:3,:),    h01_blackletter;
                  h0*e1*transpose(e_s),     H1;
                  zeros(1,3),               h1*transpose(e_s)];
%{
R1_blackletter = [eye(4), R01_blackletter_accute; zeros(3,4), R1_accute];
B1_blackletter = [h0blackletter, e_s*e1'*(-imag(ritz(3))^2); [e1;0]*e_s'*h0,B1_];
R1blackletter_inv = [eye(3), -R01blackletter*inv(R1); zeros(3), inv(R1)];
h1_blackletter = R1_blackletter*B1_blackletter*R1blackletter_inv
%}
%verify results
%Q1_blackletter * h1_blackletter
%A * Q1_blackletter(:,1:6)
%Q_ref * H_ref(:,1:6)
%A * Q_ref(:,1:6)

h1_blackletter_reduced = h1_blackletter;

for i = 1:2*s
    sc = planerot(h1_blackletter_reduced(i:i+1,i));
    G = eye(m);
    G(i,i) = sc(1,1);
    G(i,i+1) = sc(1,2);
    G(i+1,i) = sc(2,1);
    G(i+1,i+1) = sc(2,2);
    h1_blackletter_reduced = G*h1_blackletter_reduced;
    h1_blackletter_reduced(i+1:m,i) = 0;
    zeta = G*zeta;
end

y = h1_blackletter_reduced(:,1:6) \ zeta;


zeta = [beta; zeros(m-1,1)];

for i = 1:2*s
    sc = planerot(H_ref(i:i+1,i));
    G = eye(m);
    G(i,i) = sc(1,1);
    G(i,i+1) = sc(1,2);
    G(i+1,i) = sc(2,1);
    G(i+1,i+1) = sc(2,2);
    H_ref = G*H_ref;
    H_ref(i+1:m,i) = 0;
    zeta = G*zeta;
end

y_ref = H_ref(:,1:6) \ zeta;

norm(b - A*Q1_blackletter(:,1:6)*y)
norm(b - A*Q_ref(:,1:6)*y_ref)
zeta(7)

Q1_blackletter(:,1:6)*y
Q_ref(:,1:6)*y_ref
x;

% m = 9;
% x = ones(m,1);
% b = O * x;
% q1_ref = b / norm(b,2);
% 
% 
% [Q_ref, H_ref] = arnoldi(O, q1_ref, m);
% 
% 
% for i = 1:8
%     sc = planerot(H_ref(i:i+1,i));
%     G = eye(m);
%     G(i,i) = sc(1,1);
%     G(i,i+1) = sc(1,2);
%     G(i+1,i) = sc(2,1);
%     G(i+1,i+1) = sc(2,2);
%     H_ref = G*H_ref;
%     H_ref(i+1:m,i) = 0;
% end