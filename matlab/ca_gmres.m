clear all;
m = 7;
%%{
A = [0.8329    0.3181    0.6473    0.1097    0.7720    0.5254    0.5201;
    0.2564    0.1192    0.5439    0.0636    0.9329    0.5303    0.3477;
    0.6135    0.9398    0.7210    0.4046    0.9727    0.8611    0.1500;
    0.5822    0.6456    0.5225    0.4484    0.1920    0.4849    0.5861;
    0.5407    0.4795    0.9937    0.3658    0.1389    0.3935    0.2621;
    0.8699    0.6393    0.2187    0.7635    0.6963    0.6714    0.0445;
    0.2648    0.5447    0.1058    0.6279    0.0938    0.7413    0.7549];
%}
%A = rand(m);
x = ones(m,1);
b = A * x;
x0 = zeros(m,1);
q1_ref = b - A * x0;

[Q_ref, H_ref] = arnoldi(A, q1_ref, m);

eigA = eig(A);

[leja, outidx] = modified_leja(eigA, m, ones(m,1));
ritz = transpose(leja);

s = 3;
beta = norm(q1_ref);
zeta = [beta; zeros(m-1,1)];
v1 = q1_ref / beta;

B0_ = eye(3);
B0_ = [ritz(1),0,0; B0_];
B0_(2,2) = real(ritz(2));
B0_(3,3) = real(ritz(2));
B0_(2,3) = - imag(ritz(2))^2;
v2 = (A-ritz(1)*eye(m))*v1;
v3 = (A-real(ritz(2))*eye(m))*v2;
v4 = (A-real(ritz(2))*eye(m))*v3 + imag(ritz(2))^2*v2;
V0_accute = [v2, v3, v4];

V_bad = [v1, A*v1, A*A*v1, A*A*A*v1];
cond([v1, V0_accute])
cond(V_bad)

[Q0_, R0_] = qr([v1,V0_accute], 0);

R0_(1,:) = R0_(1,:) * (-1);
Q0_(:,1) = Q0_(:,1) * (-1);

R0 = R0_(1:3,1:3);
Q0 = Q0_(:, 1:3);

Q0_blackletter = Q0_;

h0_blackletter = R0_ * B0_ * inv(R0);

%verify result
A * Q0;
Q0_* h0_blackletter;
Q_ref * H_ref;

h0 = h0_blackletter(4,3);
H0 = h0_blackletter;

B1_ = eye(3);
B1_ = [ritz(4) 0 0; B1_];
B1_(2,2) = real(ritz(5));
B1_(3,3) = real(ritz(5));
B1_(2,3) = - imag(ritz(5))^2;
B1 = B1_(1:3,1:3);

v4 = Q0_(:,4);

v5 = (A - ritz(4)*eye(m))*v4;
v6 = (A - real(ritz(5))*eye(m))*v5;
v7 = (A - real(ritz(5))*eye(m))*v6 + imag(ritz(5))^2*v5;

V1_accute = [v5, v6, v7];
V_bad = [v4, A*v4, A*A*v4, A*A*A*v4];
cond([v4, V1_accute])
cond(V_bad)

R01_blackletter_accute = transpose(Q0_blackletter) * V1_accute;

V1_accute_accute = V1_accute - Q0_blackletter * R01_blackletter_accute;
%{
[Q, R] = qr(V1_accute_accute, 0);
D = diag(sign(diag(R)))
Q1_accute = Q*D
R1_accute = D*R
%}
[Q1_accute, R1_accute] = qr(V1_accute_accute, 0);

R1 = eye(3);
R1(1,:) = [1, R01_blackletter_accute(4,1), R01_blackletter_accute(4,2)];
R1(2:3,2:3) = R1_accute(1:2, 1:2);

R01_blackletter = transpose(Q0) * [v4, V1_accute(:,1:3)];
R01blackletter = transpose(Q0) * [v4, V1_accute(:,1:2)];

h01_blackletter = -h0_blackletter(1:3,1:3) * R01blackletter * inv(R1) + R01_blackletter * B1_ * inv(R1);

roh1 = R1_accute(3,3);

roh1_tilde_inv = inv(R1(3,3));

z1 = R01_blackletter_accute(4,3);
z1 = [z1;R1_accute(1:2,3)];

e_s = [0;0;1];
e1 = [1;0;0];

Q1_blackletter = [Q0_blackletter, Q1_accute];
H1 = R1 * B1 * inv(R1) + roh1_tilde_inv*1*z1*transpose(e_s) - h0 * e1 * transpose(e_s) * R01blackletter * inv(R1);

h1 = roh1_tilde_inv * roh1 * 1;

h1_blackletter = [h0_blackletter(1:3,:),    h01_blackletter;
                  h0*e1*transpose(e_s),     H1;
                  zeros(1,3),               h1*transpose(e_s)];
						
%Q1_blackletter * h1_blackletter
%A * Q1_blackletter(:,1:6)
%Q_ref * H_ref(:,1:6)
%A * Q_ref(:,1:6)

h1_blackletter_reduced = h1_blackletter;

for i = 1:m-1
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

h1_blackletter_reduced = H_ref;
zeta = [beta; zeros(m-1,1)];

for i = 1:m-1
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

y_ref = h1_blackletter_reduced(:,1:6) \ zeta;

norm(b - A*Q1_blackletter(:,1:6)*y);
norm(b - A*Q_ref(:,1:6)*y_ref);
zeta(7);

Q1_blackletter(:,1:6)*y
Q_ref(:,1:6)*y_ref;
x;

