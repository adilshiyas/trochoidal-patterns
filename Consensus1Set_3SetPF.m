% clc; clear;
% rng(0,'twister');

% Pattern generation from consensus inputs, not time-based
%2 Instances of Consensus with same parameters except phase angle

sensor_rad = 1;
holdcase = 0;
% holdcase = 1;
FigNo = 20;
PhiA = 0; PhiB = 0;
% PhiA2 = 0; PhiB2 = pi/4;
robotrad = 0.5;

k = 2; 
PatternType = 0;
CollDist = 0.5;
CommRange = 10;

%3,4,5 triplet
% B1 = (3+2*k)/(k-1);
% B3 = B1 - 3;
% B2 = 2;

% B1 = (6-k)/(k-1);
% B3 = 3 + B1;
% B2 = 2;

%5,12,13 triplet
% s1 = 3; s2 = 4; s3 = 5; %3, 4, 5
s1 = 5; s2 = 12; s3 = 13; %5, 12, 13
% s1 = 8; s2 = 15; s3 = 17; %8, 15, 17
% s1 = 7; s2 = 24; s3 = 25; %7, 24, 25
% s1 = 16; s2 = 63; s3 = 65;         %16, 63, 65
% s1 = 11; s2 = 60; s3 = 61;         %11, 60, 61

if(PatternType ==0)
    JR = k + 1;
else
    JR = k - 1;
    JR = -JR;
end

B2 = s2/2;
% B1 - B3 = s1;
B3 =  (s3 + s2 - s1 + s1*JR - JR*s2 + JR*s3)/(2*(JR-1));
B1 = B3 - s1;

scaleFactor = 1;
B1 = B1*scaleFactor;
B2 = B2*scaleFactor;
B3 = B3*scaleFactor;

B = diag([B1 B2 B3]);

a = B1/2 + B2 + B3/2;
b = (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2;
alpha3r = B3/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha3d = alpha3r;
alpha2r = (B3 - (a-b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha2d = (B3 - (a+b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1r = (B1*(a-b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1d = (B1*(a+b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));

G2 = 2*b*B1 + B1^2 - B1*B3 - 2*B2*B3;
G3 = -(2*b*B1 + B1^2 - B1*B3 -B2*B3 + 2*b*B2 + 2*B2^2 + B1*B2);
G1 = B2*B3 + 2*b*B2 + 2*B2^2 +B1*B2;

G3p = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2p = -(B1*2*b - B1^2 + B1*B3 + 2*B2*B3);
G1p = B1*B2 + 2*B2^2 + B2*B3 - B2*2*b;

G1diff = G1 - G1p;
G2diff = G2 - G2p;
G3diff = G3 - G3p;

G3Phi = B1*2*b + B2*2*b + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3;
G2Phi = -(B1*2*b + B1^2 - B1*B3 - 2*B2*B3);
G1Phi = -(B2*B1 + 2*B2^2 + B2*B3 + B2*2*b);

G3Phip = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2Phip = -(B1*2*b -B1^2 + B1*B3 + 2*B2*B3);
G1Phip = B2*B1 + 2*B2^2 + B3*B2 - B2*2*b;

D1 = [-1 0;1 -1;0 1];
% D2 = [-1 0 1;1 -1 0;0 1 -1];
% D = [1 0;-1 1;0 -1];
% D = [-1 0;1 1;0 -1];
D2 = [-1 -1;1 0;0 1];
% D3 = [-1 0;1 1;0 -1];
D3 = [-1 0;0 -1; 1 1];
% D = [1 -1]';

D4 = [-1 0 1;1 -1 0; 0 1 -1]; %Cc
% D = [1 0 1;-1 -1 0; 0 1 -1]; %Cd
% D = [1 0 -1;-1 -1 0; 0 1 1]; %Ce
% D = [1 0 1;1 1 0; 0 1 1]; %Complete undirected graph
% % D = [-1 0 -1;1 1 0; 0 -1 1]; %Cf
% D = [-1 0 1;1 -1 0;0 1 -1]; %Ca
% D = [-1 0;1 1;0 -1];

D = D1;

sigma = 1;

G = diag([0 0 0]);
% G = diag([0.1,0.5]);
% G = diag([0 0 0]);
% G = diag([3 4]);

L = D*D';
S = [0 -1; 1 0];

[vBL,eBL]=eig(B*L);
reBL = round(diag(eBL),3);
LJ = min(abs(nonzeros(reBL)));
RR = 2*pi/LJ;

% FigNo = 15;
% C01 = PolyVert(FigNo,:);

C01 = [540 225];
C01 = [1600 1000]; %k = 4
% C01 = [2400 1525];
% C01 = dscaled';
C01 = C01*(scaleFactor)^2;
% C01 = [16000 10000];
% C01 = [1000 2000];
% C01 = [16 10];
% C01 = [2000 0];
% C01 = [2000 0];
% C01 = [500,3200];
% C01 = [4300,500];
% C01 = [3800,1400];
syms x10 x20 y10 y20


eqx1 = C01(1) == sqrt((G1*x10 + G2*x20)^2 + (G1*y10 + G2*y20)^2);
eqx2 = C01(2) == sqrt((G1p*x10 + G2p*x20)^2 + (G1p*y10 + G2p*y20)^2);

eqx3 = atan2((G1Phi*y10 + G2Phi*y20),(G1Phi*x10 + G2Phi*x20)) == PhiA;
eqx4 = atan2((G1Phip*y10 + G2Phip*y20),(G1Phip*x10 + G2Phip*x20)) == PhiB;

HardSolve = solve([eqx1,eqx2,eqx3,eqx4],[x10, x20, y10, y20]);

% syms x10 x20 y10 y20
% 
% 
% eqx1 = C01(1) == sqrt((G1*x10 + G2*x20)^2 + (G1*y10 + G2*y20)^2);
% eqx2 = C01(2) == sqrt((G1p*x10 + G2p*x20)^2 + (G1p*y10 + G2p*y20)^2);
% 
% eqx3 = atan2((G1Phi*y10 + G2Phi*y20),(G1Phi*x10 + G2Phi*x20)) == PhiA2;
% eqx4 = atan2((G1Phip*y10 + G2Phip*y20),(G1Phip*x10 + G2Phip*x20)) == PhiB2;
% 
% HardSolve2 = solve([eqx1,eqx2,eqx3,eqx4],[x10, x20, y10, y20]);


% x0 = [-0.5;0.5;0.5;-0.5;1;-1];
% x0 = rand(6,1);
% x0 = [-1;0;1;0]; %2 agents consensus
% x0 = [-1;1;2;3];
% x0 = [6;8;-7;5;5;-10];
% x0 = [1; 0; -3; 0; 2; 0];
% x0 = [0.7922
%     0.9595
%     0.6557
%     0.0357
%     0.8491
%     0.9340];

%collision case
x0 = [   -0.4009
   -0.9605
    1.2003
   -0.2743
    1.6426
   -1.2726];
x0 = [  -0.5092
    0.6676
    0.3584
    0.1907
    0.3678
    0.5929];

% x0 = -5 + 10*rand(6,1);

% x0 = [1;0;4;0;10;0];

% x0 = [1;2; -3; 1; 0; 1]; %Case 1
% x0 = [1;0;-1;0;0;1]; %Case 2
% x0 = [-1;0;1;0;0;1]; %Case 3
% x0 = [-1;0;0;1;1;0]; %Case 4
% x0 = [0;1;-1;0;1;0]; %Case 5
% x0 = [0;1;1;0;-1;0]; %Case 6

% x0 = [1;0;-1;0;0;1]; %Case A
% x0 = [1;0;0;1;-1;0]; %Case B
% x0 = [-1;0;0;1;1;0]; %Case C
% x0 = [-1;0;1;0;0;1]; %Case D
% x0 = [0;1;1;0;-1;0]; %Case E
% x0 = [0;1;-1;0;1;0]; %Case F

% x0 = [1;2;4;-1;2;3];
% x0 = [1;0;-2;0;5;0];
% x0 = [4;-1;1;2;2;3];
% x0 = [2;3;1;2;4;-1];
% x0 = [1;2;2;3;4;-1];
% x0 = [2;3;4;-1;1;2];

% x0 = [0.7;1;0.5;1;-0.7;1];
% x0 = [1;0.7;1;0.5;1;-0.7];
% x0 = [0;1;-sqrt(3)/2;-0.5;sqrt(3)/2;-0.5];

% x0 = [4.4049;-1.3589;-0.7966;-2.7179;0;0]; %alpha3
% x0 = [7.9019;-1.4732;0.1787;-2.9463;0;0]; %alpha2
% x0 = [1.0723;-17.6777;-29.1054;-35.554;0;0]; %alpha1
% x0 = [-36.4276;-17.6777;-41.6053;-35.3553;0;0];
% x0 = [-36.4276;-55.1776;-41.6053;-110.3552;0;0];

% x0 = [2.1454;-35.3553;-58.2104;-70.7106;0;0];
% x0 = [2.1454;0;-58.2104;0;0;0];

l = 2;
% x0 = [-l/2;-sqrt(3)*l/4;l/2;-sqrt(3)*l/4;0;sqrt(3)*l/4];

L12 = 1;
L23 = L12;
The13 = 0*pi/3;
The23 = 6*pi/6;
% x0 = [1;0;2;0;0;0];
% x0 = [L12*cos(The13);L12*sin(The13);L23*cos(The23);L23*sin(The23);0;0];

% x0 = [-28.1680;0;28.1680;0;0;0];
x0 = [26.7963;8.6830;10.7185;26.0490;0;0];
x0 = [2.4772;0.5031;0.9909;1.5094;0;0];
x0 = [-2.4772;0.5031;-0.9909;1.5094;0;0];
x0 = [-91/36;0; 65/6;0;0;0];
x0 = [-23;0;-35;0;0;0];
x0 = [13;13;5.2;39;0;0];
% x0 = [1;0;3;0;0;0]
x0 = [-7.5;0;-2.5;0;0;0];
x0 = [-1.09;0;-0.363;0;0;0];
x0 = [-20;0;-20;0;0;0];
x0 = [-2;0;-14;0;0;0];
x0 = [-13;0;-11;0;0;0];
x0 = [-16;0;-7;0;0;0];
x0 = [-11;0;-12;0;0;0];
x0 = [-17.714;0;-9;0;0;0];
x0 = [-11.5;0;-10.5;0;0;0];
x0 = [-8.735;0;-11.834;0;0;0];
x0 = [-3.14;0;-1.04;0;0;0];
x0 = [-0.0191;0;-0.065207;0;0;0];
x0 = [1;0;-3;0;0;0];

x0 = [round(HardSolve.x10,4);round(HardSolve.y10,4);round(HardSolve.x20,4);round(HardSolve.y20,4);0;0];
% x0_2 = [round(HardSolve2.x10,4);round(HardSolve2.y10,4);round(HardSolve2.x20,4);round(HardSolve2.y20,4);0;0];

PhiD_Set2 = pi/2;
PhiD_Set3 = pi;
PhiD_Set4 = 3*pi/2;

% x0_2 = ParametricFormConsensus2(vBL,eBL,PhiD_Set1/LJ,S,x0);

% x0 = [1;0;-3;0;0;0];
% x0 = [1;0;2.5;0;0;0];
% x0 = [1;0;2;0;0;0];

x0_Orig = x0;
% x01 = -1.5 + 3*rand(2,1);
% x02 = -1.5 + 3*rand(2,1);
% while(1)
%     if(norm([x01(1) x01(2)] - [x02(1) x02(2)])<=0.5)
%         x02 = -1.5 + 3*rand(2,1);
%     else
%         break
%     end
% end
% 
% x03 = -1.5 + 3*rand(2,1);

% while(1)
%     if(norm([x03(1) x03(2)] - [x02(1) x02(2)])<=0.5)||(norm([x01(1) x01(2)] - [x03(1) x03(2)])<=0.5)
%         x03 = -1.5 + 3*rand(2,1);
%     else
%         break
%     end
% end
% x0 = [x01;x02;x03];

M = -(kron((G*L),eye(2,2))-kron((B*L),S));
V=null(M');
M1 = [V(:,1)'*(kron(ones(3,1),eye(2,2))); V(:,2)'*(kron(ones(3,1),eye(2,2)))];
M2 = [V(:,1)'*x0; V(:,2)'*x0];
xInf = M1\M2;

x0 = x0 - repmat(xInf,3,1);

x0_2 = ParametricFormConsensus2(vBL,eBL,PhiD_Set2/LJ,S,x0);
x0_3 = ParametricFormConsensus2(vBL,eBL,PhiD_Set3/LJ,S,x0);
x0_4 = ParametricFormConsensus2(vBL,eBL,PhiD_Set4/LJ,S,x0);

% x0_2 = x0_2 - repmat(xInf_2,3,1);

M2 = [V(:,1)'*x0; V(:,2)'*x0];
xInfR_1 = M1\M2;

r1 = alpha1r*C01(1)/(k+1);
d1 = alpha1d*C01(2);
r2 = alpha2r*C01(1)/(k+1);
d2 = alpha2d*C01(2);
r3 = alpha3r*C01(1)/(k+1);
d3 = alpha3d*C01(2);


dT = 0.001;
tend = RR;
t = 0:dT:tend; 

dist12 = zeros(1,length(t));
dist13 = dist12;
% dist14 = dist12;
dist15 = dist12;
dist16 = dist12;
% dist17 = dist12;
dist18 = dist12;
dist19 = dist12;
% dist1_10 = dist12;
dist1_11 = dist12;
dist1_12 = dist12;
dist23 = dist12;
dist24 = dist12;
% dist25 = dist12;
dist26 = dist12;
dist27 = dist12;
% dist28 = dist12;
dist29 = dist12;
dist2_10 = dist12;
% dist2_11 = dist12;
dist2_12 = dist12;
dist34 = dist12;
dist35 = dist12;
% dist36 = dist12;
dist37 = dist12;
dist38 = dist12;
% dist39 = dist12;
dist3_10 = dist12;
dist3_11 = dist12;
% dist3_12 = dist12;
% dist45 = dist12;
% dist46 = dist12;
% dist47 = dist12;
dist48 = dist12;
dist49 = dist12;
% dist4_10 = dist12;
dist4_11 = dist12;
dist4_12 = dist12;
% dist56 = dist12;
dist57 = dist12;
% dist58 = dist12;
dist59 = dist12;
dist5_10 = dist12;
% dist5_11 = dist12;
dist5_12 = dist12;
dist67 = dist12;
dist68 = dist12;
% dist69 = dist12;
dist6_10 = dist12;
dist6_11 = dist12;
% dist6_12 = dist12;
% dist78 = dist12;
% dist79 = dist12;
% dist7_10 = dist12;
dist7_11 = dist12;
dist7_12 = dist12;
% dist89 = dist12;
dist8_10 = dist12;
% dist8_11 = dist12;
dist8_12 = dist12;
dist9_10 = dist12;
dist9_11 = dist12;
% dist9_12 = dist12;
% dist10_11 = dist12;
% dist10_12 = dist12;



dist1 = zeros(1,length(t));
dist3 = dist1;
dist2 = dist1;

x = zeros(6,length(t)); u = x;

x_2 = x; u_2 = x;

i = 1;
x(:,i) = x0;
z = kron(D',eye(2))*x(:,i);
u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

x_2(:,i) = x0_2;
z_2 = kron(D',eye(2))*x_2(:,i);
u_2(:,i) = -(kron(G*D,eye(2,2)))*z_2 + kron((B*D),S)*z_2;

x_3(:,i) = x0_3;
z_3 = kron(D',eye(2))*x_3(:,i);
u_3(:,i) = -(kron(G*D,eye(2,2)))*z_3 + kron((B*D),S)*z_3;

x_4(:,i) = x0_4;
z_4 = kron(D',eye(2))*x_4(:,i);
u_4(:,i) = -(kron(G*D,eye(2,2)))*z_4 + kron((B*D),S)*z_4;


% dist12(i) = norm([x(1,i) x(2,i)] - [x(3,i) x(4,i)]);
% dist31(i) = norm([x(1,i) x(2,i)] - [x(5,i) x(6,i)]);
% dist23(i) = norm([x(5,i) x(6,i)] - [x(3,i) x(4,i)]);

dist12(i) = norm([x(1,i) x(2,i)] - [x(3,i) x(4,i)]);
dist13(i) = norm([x(1,i) x(2,i)] - [x(5,i) x(6,i)]);
% dist14(i) = norm([x(1,i) x(2,i)] - [x_2(1,i) x_2(2,i)]);
dist15(i) = norm([x(1,i) x(2,i)] - [x_2(3,i) x_2(4,i)]);
dist16(i) = norm([x(1,i) x(2,i)] - [x_2(5,i) x_2(6,i)]);
% dist17(i) = norm([x(1,i) x(2,i)] - [x_3(1,i) x_3(2,i)]);
dist18(i) = norm([x(1,i) x(2,i)] - [x_3(3,i) x_3(4,i)]);
dist19(i) = norm([x(1,i) x(2,i)] - [x_3(5,i) x_3(6,i)]);
% dist1_10(i) = norm([x(1,i) x(2,i)] - [x_4(1,i) x_4(2,i)]);
dist1_11(i) = norm([x(1,i) x(2,i)] - [x_4(3,i) x_4(4,i)]);
dist1_12(i) = norm([x(1,i) x(2,i)] - [x_4(5,i) x_4(6,i)]);

dist23(i) = norm([x(5,i) x(6,i)] - [x(3,i) x(4,i)]);
dist24(i) = norm([x(3,i) x(4,i)] - [x_2(1,i) x_2(2,i)]);
% dist25(i) = norm([x(3,i) x(4,i)] - [x_2(3,i) x_2(4,i)]);
dist26(i) = norm([x(3,i) x(4,i)] - [x_2(5,i) x_2(6,i)]);
dist27(i) = norm([x(3,i) x(4,i)] - [x_3(1,i) x_3(2,i)]);
% dist28(i) = norm([x(3,i) x(4,i)] - [x_3(3,i) x_3(4,i)]);
dist29(i) = norm([x(3,i) x(4,i)] - [x_3(5,i) x_3(6,i)]);
dist2_10(i) = norm([x(3,i) x(4,i)] - [x_4(1,i) x_4(2,i)]);
% dist2_11(i) = norm([x(3,i) x(4,i)] - [x_4(3,i) x_4(4,i)]);
dist2_12(i) = norm([x(3,i) x(4,i)] - [x_4(5,i) x_4(6,i)]);

dist34(i) = norm([x(5,i) x(6,i)] - [x_2(1,i) x_2(2,i)]);
dist35(i) = norm([x(5,i) x(6,i)] - [x_2(3,i) x_2(4,i)]);
% dist36(i) = norm([x(5,i) x(6,i)] - [x_2(5,i) x_2(6,i)]);
dist37(i) = norm([x(5,i) x(6,i)] - [x_3(1,i) x_3(2,i)]);
dist38(i) = norm([x(5,i) x(6,i)] - [x_3(3,i) x_3(4,i)]);
% dist39(i) = norm([x(5,i) x(6,i)] - [x_3(5,i) x_3(6,i)]);
dist3_10(i) = norm([x(5,i) x(6,i)] - [x_4(1,i) x_4(2,i)]);
dist3_11(i) = norm([x(5,i) x(6,i)] - [x_4(3,i) x_4(4,i)]);
% dist3_12(i) = norm([x(5,i) x(6,i)] - [x_4(5,i) x_4(6,i)]);

% dist45(i) = norm([x_2(1,i) x_2(2,i)] - [x_2(3,i) x_2(4,i)]);
% dist46(i) = norm([x_2(1,i) x_2(2,i)] - [x_2(5,i) x_2(6,i)]);
% dist47(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(1,i) x_3(2,i)]);
dist48(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(3,i) x_3(4,i)]);
dist49(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(5,i) x_3(6,i)]);
% dist4_10(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(1,i) x_4(2,i)]);
dist4_11(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(3,i) x_4(4,i)]);
dist4_12(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(5,i) x_4(6,i)]);

% dist56(i) = norm([x_2(3,i) x_2(4,i)] - [x_2(5,i) x_2(6,i)]);
dist57(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(1,i) x_3(2,i)]);
% dist58(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(3,i) x_3(4,i)]);
dist59(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(5,i) x_3(6,i)]);
dist5_10(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(1,i) x_4(2,i)]);
% dist5_11(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(3,i) x_4(4,i)]);
dist5_12(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(5,i) x_4(6,i)]);

dist67(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(1,i) x_3(2,i)]);
dist68(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(3,i) x_3(4,i)]);
% dist69(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(5,i) x_3(6,i)]);
dist6_10(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(1,i) x_4(2,i)]);
dist6_11(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(3,i) x_4(4,i)]);
% dist6_12(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(5,i) x_4(6,i)]);

% dist78(i) = norm([x_3(1,i) x_3(2,i)] - [x_3(3,i) x_3(4,i)]);
% dist79(i) = norm([x_3(1,i) x_3(2,i)] - [x_3(5,i) x_3(6,i)]);
% dist7_10(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(1,i) x_4(2,i)]);
dist7_11(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(3,i) x_4(4,i)]);
dist7_12(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(5,i) x_4(6,i)]);

% dist89(i) = norm([x_3(3,i) x_3(4,i)] - [x_3(5,i) x_3(6,i)]);
dist8_10(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(1,i) x_4(2,i)]);
% dist8_11(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(3,i) x_4(4,i)]);
dist8_12(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(5,i) x_4(6,i)]);

dist9_10(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(1,i) x_4(2,i)]);
dist9_11(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(3,i) x_4(4,i)]);
% dist9_12(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(5,i) x_4(6,i)]);

% dist10_11(i) = norm([x_4(1,i) x_4(2,i)] - [x_4(3,i) x_4(4,i)]);
% dist10_12(i) = norm([x_4(1,i) x_4(2,i)] - [x_4(5,i) x_4(6,i)]);

% dist11_12(i) = norm([x_4(3,i) x_4(4,i)] - [x_4(5,i) x_4(6,i)]);


coll_flag = 0;
dist_thresh = 1;
retval = 1;
switchiter = [];



for i = 2:length(t)    


    x(:,i) = x(:,i-1) + dT*eye(6,6)*(u(:,i-1));
    x_2(:,i) = x_2(:,i-1) + dT*eye(6,6)*(u_2(:,i-1));
    x_3(:,i) = x_3(:,i-1) + dT*eye(6,6)*(u_3(:,i-1));
    x_4(:,i) = x_4(:,i-1) + dT*eye(6,6)*(u_4(:,i-1));
   

    z = kron(D',eye(2))*x(:,i);
    u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

    z_2 = kron(D',eye(2))*x_2(:,i);
    u_2(:,i) = -(kron(G*D,eye(2,2)))*z_2 + kron((B*D),S)*z_2;

    z_3 = kron(D',eye(2))*x_3(:,i);
    u_3(:,i) = -(kron(G*D,eye(2,2)))*z_3 + kron((B*D),S)*z_3;

    z_4 = kron(D',eye(2))*x_4(:,i);
    u_4(:,i) = -(kron(G*D,eye(2,2)))*z_4 + kron((B*D),S)*z_4;
    
    dist12(i) = norm([x(1,i) x(2,i)] - [x(3,i) x(4,i)]);
    dist13(i) = norm([x(1,i) x(2,i)] - [x(5,i) x(6,i)]);
%     dist14(i) = norm([x(1,i) x(2,i)] - [x_2(1,i) x_2(2,i)]);
    dist15(i) = norm([x(1,i) x(2,i)] - [x_2(3,i) x_2(4,i)]);
    dist16(i) = norm([x(1,i) x(2,i)] - [x_2(5,i) x_2(6,i)]);
%     dist17(i) = norm([x(1,i) x(2,i)] - [x_3(1,i) x_3(2,i)]);
    dist18(i) = norm([x(1,i) x(2,i)] - [x_3(3,i) x_3(4,i)]);
    dist19(i) = norm([x(1,i) x(2,i)] - [x_3(5,i) x_3(6,i)]);
%     dist1_10(i) = norm([x(1,i) x(2,i)] - [x_4(1,i) x_4(2,i)]);
    dist1_11(i) = norm([x(1,i) x(2,i)] - [x_4(3,i) x_4(4,i)]);
    dist1_12(i) = norm([x(1,i) x(2,i)] - [x_4(5,i) x_4(6,i)]);
    
    dist23(i) = norm([x(5,i) x(6,i)] - [x(3,i) x(4,i)]);
    dist24(i) = norm([x(3,i) x(4,i)] - [x_2(1,i) x_2(2,i)]);
%     dist25(i) = norm([x(3,i) x(4,i)] - [x_2(3,i) x_2(4,i)]);
    dist26(i) = norm([x(3,i) x(4,i)] - [x_2(5,i) x_2(6,i)]);
    dist27(i) = norm([x(3,i) x(4,i)] - [x_3(1,i) x_3(2,i)]);
%     dist28(i) = norm([x(3,i) x(4,i)] - [x_3(3,i) x_3(4,i)]);
    dist29(i) = norm([x(3,i) x(4,i)] - [x_3(5,i) x_3(6,i)]);
    dist2_10(i) = norm([x(3,i) x(4,i)] - [x_4(1,i) x_4(2,i)]);
%     dist2_11(i) = norm([x(3,i) x(4,i)] - [x_4(3,i) x_4(4,i)]);
    dist2_12(i) = norm([x(3,i) x(4,i)] - [x_4(5,i) x_4(6,i)]);
    
    dist34(i) = norm([x(5,i) x(6,i)] - [x_2(1,i) x_2(2,i)]);
    dist35(i) = norm([x(5,i) x(6,i)] - [x_2(3,i) x_2(4,i)]);
%     dist36(i) = norm([x(5,i) x(6,i)] - [x_2(5,i) x_2(6,i)]);
    dist37(i) = norm([x(5,i) x(6,i)] - [x_3(1,i) x_3(2,i)]);
    dist38(i) = norm([x(5,i) x(6,i)] - [x_3(3,i) x_3(4,i)]);
%     dist39(i) = norm([x(5,i) x(6,i)] - [x_3(5,i) x_3(6,i)]);
    dist3_10(i) = norm([x(5,i) x(6,i)] - [x_4(1,i) x_4(2,i)]);
    dist3_11(i) = norm([x(5,i) x(6,i)] - [x_4(3,i) x_4(4,i)]);
%     dist3_12(i) = norm([x(5,i) x(6,i)] - [x_4(5,i) x_4(6,i)]);
    
%     dist45(i) = norm([x_2(1,i) x_2(2,i)] - [x_2(3,i) x_2(4,i)]);
%     dist46(i) = norm([x_2(1,i) x_2(2,i)] - [x_2(5,i) x_2(6,i)]);
%     dist47(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(1,i) x_3(2,i)]);
    dist48(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(3,i) x_3(4,i)]);
    dist49(i) = norm([x_2(1,i) x_2(2,i)] - [x_3(5,i) x_3(6,i)]);
%     dist4_10(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(1,i) x_4(2,i)]);
    dist4_11(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(3,i) x_4(4,i)]);
    dist4_12(i) = norm([x_2(1,i) x_2(2,i)] - [x_4(5,i) x_4(6,i)]);
    
%     dist56(i) = norm([x_2(3,i) x_2(4,i)] - [x_2(5,i) x_2(6,i)]);
    dist57(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(1,i) x_3(2,i)]);
%     dist58(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(3,i) x_3(4,i)]);
    dist59(i) = norm([x_2(3,i) x_2(4,i)] - [x_3(5,i) x_3(6,i)]);
    dist5_10(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(1,i) x_4(2,i)]);
%     dist5_11(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(3,i) x_4(4,i)]);
    dist5_12(i) = norm([x_2(3,i) x_2(4,i)] - [x_4(5,i) x_4(6,i)]);
    
    dist67(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(1,i) x_3(2,i)]);
    dist68(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(3,i) x_3(4,i)]);
%     dist69(i) = norm([x_2(5,i) x_2(6,i)] - [x_3(5,i) x_3(6,i)]);
    dist6_10(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(1,i) x_4(2,i)]);
    dist6_11(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(3,i) x_4(4,i)]);
%     dist6_12(i) = norm([x_2(5,i) x_2(6,i)] - [x_4(5,i) x_4(6,i)]);
    
%     dist78(i) = norm([x_3(1,i) x_3(2,i)] - [x_3(3,i) x_3(4,i)]);
%     dist79(i) = norm([x_3(1,i) x_3(2,i)] - [x_3(5,i) x_3(6,i)]);
%     dist7_10(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(1,i) x_4(2,i)]);
    dist7_11(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(3,i) x_4(4,i)]);
    dist7_12(i) = norm([x_3(1,i) x_3(2,i)] - [x_4(5,i) x_4(6,i)]);
    
%     dist89(i) = norm([x_3(3,i) x_3(4,i)] - [x_3(5,i) x_3(6,i)]);
    dist8_10(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(1,i) x_4(2,i)]);
%     dist8_11(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(3,i) x_4(4,i)]);
    dist8_12(i) = norm([x_3(3,i) x_3(4,i)] - [x_4(5,i) x_4(6,i)]);
    
    dist9_10(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(1,i) x_4(2,i)]);
    dist9_11(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(3,i) x_4(4,i)]);
%     dist9_12(i) = norm([x_3(5,i) x_3(6,i)] - [x_4(5,i) x_4(6,i)]);
    
%     dist10_11(i) = norm([x_4(1,i) x_4(2,i)] - [x_4(3,i) x_4(4,i)]);
%     dist10_12(i) = norm([x_4(1,i) x_4(2,i)] - [x_4(5,i) x_4(6,i)]);
%     
%     dist11_12(i) = norm([x_4(3,i) x_4(4,i)] - [x_4(5,i) x_4(6,i)]);
    
end


%Finding angle b/w x0 wrt origin

% dot13 = dot([x0(1) x0(2)],[x0(5) x0(6)]);
% mag1 = norm([x0(1) x0(2)]);
% mag3 = norm([x0(5) x0(6)]);
% ang13 = acos(dot13/(mag1*mag3));
% 
% LOS1unwrap = unwrap(LOS);
% 
% I1unwrap = zeros(1,length(t));
% for m = 2:(length(t))
%     I1unwrap(m) = ((LOS1unwrap(m) - LOS1unwrap(m-1))/dT)/((sigma^2)*(R1(m)^2));
% 
% end
% 
% 
% I1unwrap_tot = 0;
% for i = 1:length(I1)
%     I1unwrap_tot = I1unwrap_tot + I1unwrap(i);
% end




% I = I/(sigma)^2;
% x = x(:,1:i);

Acl = kron((D'*G*D),eye(2)) + kron((D'*B*D),S);
rootbl = sqrt(B(1,1)^2 - 2*B(1,1)*B(3,3) + 4*(B(2,2)^2) + B(3,3)^2)/2;
nonrootbl = B(1,1)/2 + B(2,2) + B(3,3)/2;

eigenval = eig(B*L);
eir = [];
for ei = 1:length(eigenval)
    if(abs(eigenval(ei))>0.1)
        eir = [eir eigenval(ei)];
    end
end

% eirm = [abs(eir(1)) abs(eir(2))]; 
% eiratio = max(eirm)/min(eirm);
% 
% if sign(eir(1)/eir(2))==-1
%     k = eiratio + 1;
% else
%     k = eiratio - 1;
% end

x0x = [x0(1) x0(3) x0(5)];
x0y = [x0(2) x0(4) x0(6)];
x10 = x0x(1); x20 = x0x(2); x30 = x0x(3);
y10 = x0y(1); y20 = x0y(2); y30 = x0y(3);

B1 = B(1,1);
B2 = B(2,2);
B3 = B(3,3);

BL = B*L;
[vBL,eBL] = eig(BL);
vBLinv = inv(vBL);

ttest = 0.001;

c3r = (((B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2);
c1r = -((B1*B2 + B1*B3)/(B2*B3) - (B1*(B1/2 + B2 + B3/2 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2))/(B2*B3))*(((B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2);
c2r = -((B1/2 + B2 + B3/2 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2)/B3 - 1)*(((B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2);


c3r2 = (((B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 - B1*B3 - 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2);

c1d = -(((B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2)*((B1*B2 + B1*B3)/(B2*B3) - (B1*(B1/2 + B2 + B3/2 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2))/(B2*B3));
c2d = -(((B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2)*((B1/2 + B2 + B3/2 + (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2)/B3 - 1);
c3d = (((B3*conj(x30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(x20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(x10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2 + ((B3*conj(y30)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) - (B3*conj(y20)*(B1*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) - B1^2 + B1*B3 + 2*B2*B3))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))) + (B2*B3*conj(y10)*(B1 + 2*B2 + B3 - (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)))/(2*(B1*B2*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B1*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2) + B2*B3*(B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2))))^2)^(1/2);

% c1d =

c11 = vBL(1,1)*sqrt( (vBLinv(1,:)*x0y')^2 + (vBLinv(1,:)*x0x')^2);
c12 = vBL(1,2)*sqrt( (vBLinv(2,:)*x0y')^2 + (vBLinv(2,:)*x0x')^2);
c13 = vBL(1,3)*sqrt( (vBLinv(3,:)*x0y')^2 + (vBLinv(3,:)*x0x')^2);
psix11 = atan2(vBLinv(1,:)*x0y',vBLinv(1,:)*x0x');
psix12 = atan2(vBLinv(2,:)*x0y',vBLinv(2,:)*x0x');
psix13 = atan2(vBLinv(3,:)*x0y',vBLinv(3,:)*x0x');

x1 = c11*cos(eBL(1,1)*ttest + psix11) + c12*cos(eBL(2,2)*ttest + psix12) + c13*cos(eBL(3,3)*ttest + psix13);
psiy11 = atan2(-vBLinv(1,:)*x0x',vBLinv(1,:)*x0y');
psiy12 = atan2(-vBLinv(2,:)*x0x',vBLinv(2,:)*x0y');
psiy13 = atan2(-vBLinv(3,:)*x0x',vBLinv(3,:)*x0y');
y1 = c11*cos(eBL(1,1)*ttest + psiy11) + c12*cos(eBL(2,2)*ttest + psiy12) + c13*cos(eBL(3,3)*ttest + psiy13);

c21 = vBL(2,1)*sqrt( (vBLinv(1,:)*x0y')^2 + (vBLinv(1,:)*x0x')^2);
c22 = vBL(2,2)*sqrt( (vBLinv(2,:)*x0y')^2 + (vBLinv(2,:)*x0x')^2);
c23 = vBL(2,3)*sqrt( (vBLinv(3,:)*x0y')^2 + (vBLinv(3,:)*x0x')^2);
psix21 = atan2(vBLinv(1,:)*x0y',vBLinv(1,:)*x0x');
psix22 = atan2(vBLinv(2,:)*x0y',vBLinv(2,:)*x0x');
psix23 = atan2(vBLinv(3,:)*x0y',vBLinv(3,:)*x0x');

x2 = c21*cos(eBL(1,1)*ttest + psix21) + c22*cos(eBL(2,2)*ttest + psix22) + c23*cos(eBL(3,3)*ttest + psix23);
psiy21 = atan2(-vBLinv(1,:)*x0x',vBLinv(1,:)*x0y');
psiy22 = atan2(-vBLinv(2,:)*x0x',vBLinv(2,:)*x0y');
psiy23 = atan2(-vBLinv(3,:)*x0x',vBLinv(3,:)*x0y');
y2 = c21*cos(eBL(1,1)*ttest + psiy21) + c22*cos(eBL(2,2)*ttest + psiy22) + c23*cos(eBL(3,3)*ttest + psiy23);

c31 = vBL(3,1)*sqrt( (vBLinv(1,:)*x0y')^2 + (vBLinv(1,:)*x0x')^2);
c32 = vBL(3,2)*sqrt( (vBLinv(2,:)*x0y')^2 + (vBLinv(2,:)*x0x')^2);
c33 = vBL(3,3)*sqrt( (vBLinv(3,:)*x0y')^2 + (vBLinv(3,:)*x0x')^2);
psix31 = atan2(vBLinv(1,:)*x0y',vBLinv(1,:)*x0x');
psix32 = atan2(vBLinv(2,:)*x0y',vBLinv(2,:)*x0x');
psix33 = atan2(vBLinv(3,:)*x0y',vBLinv(3,:)*x0x');

x3 = c31*cos(eBL(1,1)*ttest + psix31) + c32*cos(eBL(2,2)*ttest + psix32) + c33*cos(eBL(3,3)*ttest + psix33);
psiy31 = atan2(-vBLinv(1,:)*x0x',vBLinv(1,:)*x0y');
psiy32 = atan2(-vBLinv(2,:)*x0x',vBLinv(2,:)*x0y');
psiy33 = atan2(-vBLinv(3,:)*x0x',vBLinv(3,:)*x0y');
y3 = c31*cos(eBL(1,1)*ttest + psiy31) + c32*cos(eBL(2,2)*ttest + psiy32) + c33*cos(eBL(3,3)*ttest + psiy33);

a1dist = abs(abs(c1r) - abs(c1d));
a2dist = abs(abs(c2r) - abs(c2d));
a3dist = abs(abs(c3r) - abs(c3d));

disp([B1 B2 B3]);
disp([a1dist a2dist a3dist]);

a = B1/2 + B2 + B3/2;
b = (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2;
alpha3r = B3/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha3d = alpha3r;
alpha2r = (B3 - (a-b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha2d = (B3 - (a+b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1r = (B1*(a-b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1d = (B1*(a+b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));

G2 = 2*b*B1 + B1^2 - B1*B3 - 2*B2*B3;
G3 = -(2*b*B1 + B1^2 - B1*B3 -B2*B3 + 2*b*B2 + 2*B2^2 + B1*B2);
G1 = B2*B3 + 2*b*B2 + 2*B2^2 +B1*B2;

G3p = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2p = -(B1*2*b - B1^2 + B1*B3 + 2*B2*B3);
G1p = B1*B2 + 2*B2^2 + B2*B3 - B2*2*b;

G1diff = G1 - G1p;
G2diff = G2 - G2p;
G3diff = G3 - G3p;

G3Phi = B1*2*b + B2*2*b + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3;
G2Phi = -(B1*2*b + B1^2 - B1*B3 - 2*B2*B3);
G1Phi = -(B2*B1 + 2*B2^2 + B2*B3 + B2*2*b);

G3Phip = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2Phip = -(B1*2*b -B1^2 + B1*B3 + 2*B2*B3);
G1Phip = B2*B1 + 2*B2^2 + B3*B2 - B2*2*b;


cres = [c13 c11;c23 c21;c33 c31];
% FigNo = 1;
figure(FigNo*10); clf(FigNo*10);
plot(x(1,:),x(2,:))
hold on
plot(x_2(1,:),x_2(2,:),'c')

plot(x(3,:),x(4,:),'r')
plot(x_2(3,:),x_2(4,:),'m')

plot(x(5,:),x(6,:),'g')
plot(x_2(5,:),x_2(6,:),'k')



plot(x(1,1),x(2,1),'ob')
plot(x(3,1),x(4,1),'or')
plot(x(5,1),x(6,1),'og')

plot(x_2(1,1),x_2(2,1),'oc')
plot(x_2(3,1),x_2(4,1),'om')
plot(x_2(5,1),x_2(6,1),'ok')
axis equal
% plot(xInf(1),xInf(2),'ok')
% plot(xInf2(1),xInf2(2),'ok')
Dg = 1;
thguard = 0:0.01:2*pi;
plot(Dg*cos(thguard(1:end)),Dg*sin(thguard(1:end)),'k')




% figure(5); clf(5);
% subplot(3,1,1)
% plot(t,dist12)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% subplot(3,1,2)
% plot(t,dist23)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% subplot(3,1,3)
% plot(t,dist13)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')

tcirc = 0:0.1:2*pi;
tot_area = 0;

figure(6); clf(6);
subplot(3,1,1)
plot(t,dist1)
subplot(3,1,2)
plot(t,dist2)
subplot(3,1,3)
plot(t,dist3)

ulv = zeros(3,length(u));

for i = 1:length(u)
    ulv(1,i) = sqrt(u(1,i)^2 + u(2,i)^2);
    ulv(2,i) = sqrt(u(3,i)^2 + u(4,i)^2);
    ulv(3,i) = sqrt(u(5,i)^2 + u(6,i)^2);
end

figure(7); clf(7);
subplot(3,1,1)
plot(t,ulv(1,:))
subplot(3,1,2)
plot(t,ulv(2,:))
subplot(3,1,3)
plot(t,ulv(3,:))

% figure(8); clf(8);
% subplot(5,3,1)
% plot(t,dist12)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist12')
% subplot(5,3,2)
% plot(t,dist13)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist13')
% subplot(5,3,3)
% plot(t,dist14)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist14')
% subplot(5,3,4)
% plot(t,dist15)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist15')
% subplot(5,3,5)
% plot(t,dist16)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist16')
% subplot(5,3,6)
% plot(t,dist23)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist23')
% subplot(5,3,7)
% plot(t,dist24)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist24')
% subplot(5,3,8)
% plot(t,dist25)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist25')
% subplot(5,3,9)
% plot(t,dist26)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist26')
% subplot(5,3,10)
% plot(t,dist34)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist34')
% subplot(5,3,11)
% plot(t,dist35)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist35')
% subplot(5,3,12)
% plot(t,dist36)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist36')
% subplot(5,3,13)
% plot(t,dist45)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist45')
% subplot(5,3,14)
% plot(t,dist46)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist46')
% subplot(5,3,15)
% plot(t,dist56)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist56')
% 
% figure(9); clf(9);
% subplot(3,3,1)
% plot(t,dist14)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist14')
% subplot(3,3,2)
% plot(t,dist15)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist15')
% subplot(3,3,3)
% plot(t,dist16)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist16')
% 
% subplot(3,3,4)
% plot(t,dist24)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist24')
% subplot(3,3,5)
% plot(t,dist25)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist25')
% subplot(3,3,6)
% plot(t,dist26)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist26')
% subplot(3,3,7)
% plot(t,dist34)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist34')
% subplot(3,3,8)
% plot(t,dist35)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist35')
% subplot(3,3,9)
% plot(t,dist36)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist36')
% 
% figure(9); clf(9);
% subplot(3,3,1)
% plot(t,dist14)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist14')
% subplot(3,3,2)
% plot(t,dist15)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist15')
% subplot(3,3,3)
% plot(t,dist16)
% hold on
% plot([0,tend],[CollDist,CollDist],'r')
% plot([0,tend],[CommRange,CommRange],'g')
% title('dist16')

figure(100); clf(100);
subplot(9,4,1)
plot(t,dist15)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist15')
subplot(9,4,2)
plot(t,dist16)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist16')
subplot(9,4,3)
plot(t,dist18)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist18')
subplot(9,4,4)
plot(t,dist19)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist19')
subplot(9,4,5)
plot(t,dist1_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist1_11')
subplot(9,4,6)
plot(t,dist1_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist1_12')

subplot(9,4,7)
plot(t,dist24)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist24')
subplot(9,4,8)
plot(t,dist26)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist26')
subplot(9,4,9)
plot(t,dist27)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist27')
subplot(9,4,10)
plot(t,dist29)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist29')
subplot(9,4,11)
plot(t,dist2_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist2_10')
subplot(9,4,12)
plot(t,dist2_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist2_12')

subplot(9,4,13)
plot(t,dist34)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist34')
subplot(9,4,14)
plot(t,dist35)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist35')
subplot(9,4,15)
plot(t,dist37)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist37')
subplot(9,4,16)
plot(t,dist38)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist38')
subplot(9,4,17)
plot(t,dist3_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist3_10')
subplot(9,4,18)
plot(t,dist3_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist3_11')

subplot(9,4,19)
plot(t,dist48)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist48')
subplot(9,4,20)
plot(t,dist49)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist49')
subplot(9,4,21)
plot(t,dist4_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist4_11')
subplot(9,4,22)
plot(t,dist4_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist4_12')

subplot(9,4,23)
plot(t,dist57)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist57')
subplot(9,4,24)
plot(t,dist59)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist59')
subplot(9,4,25)
plot(t,dist5_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist5_10')
subplot(9,4,26)
plot(t,dist5_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist5_12')

subplot(9,4,27)
plot(t,dist67)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist67')
subplot(9,4,28)
plot(t,dist68)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist68')
subplot(9,4,29)
plot(t,dist6_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist6_10')
subplot(9,4,30)
plot(t,dist6_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist6_11')

subplot(9,4,31)
plot(t,dist7_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist7_11')
subplot(9,4,32)
plot(t,dist7_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist7_12')

subplot(9,4,33)
plot(t,dist8_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist8_10')
subplot(9,4,34)
plot(t,dist8_12)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist8_12')

subplot(9,4,35)
plot(t,dist9_10)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist9_10')
subplot(9,4,36)
plot(t,dist9_11)
hold on
plot([0,tend],[CollDist,CollDist],'r')
plot([0,tend],[CommRange,CommRange],'g')
title('dist9_11')
% ar1 = polyshape(sensor_rad*cos(tcirc)+x(1,1),sensor_rad*sin(tcirc)+x(2,1));
% ar2 = polyshape(sensor_rad*cos(tcirc)+x(3,1),sensor_rad*sin(tcirc)+x(4,1));
% ar3 = polyshape(sensor_rad*cos(tcirc)+x(5,1),sensor_rad*sin(tcirc)+x(6,1));
% 
% ar1tot = polyshape(0,0);
% ar2tot = polyshape(0,0);
% ar3tot = polyshape(0,0);
% 
% for i = 2:10:length(t)
% 
%     ar1 = polyshape(sensor_rad*cos(tcirc)+x(1,i),sensor_rad*sin(tcirc)+x(2,i));
%     ar2 = polyshape(sensor_rad*cos(tcirc)+x(3,i),sensor_rad*sin(tcirc)+x(4,i));
%     ar3 = polyshape(sensor_rad*cos(tcirc)+x(5,i),sensor_rad*sin(tcirc)+x(6,i));
% 
%     ar1tot = union(ar1,ar1tot);
%     ar2tot = union(ar2,ar2tot);
%     ar3tot = union(ar3,ar3tot);
% 
% end
% 
% 
% guardcirc = polyshape(Dg*cos(thguard(1:end)),Dg*sin(thguard(1:end)));
% guardarea = area(guardcirc);
% ar1totg = subtract(ar1tot,intersect(ar1tot,guardcirc));
% ar2totg = subtract(ar2tot,intersect(ar2tot,guardcirc));
% ar3totg = subtract(ar3tot,intersect(ar3tot,guardcirc));
% 
% artot = union(ar1totg,ar2totg);
% artot = union(artot,ar3totg);
% totarea = area(artot);
% 
% % figure(FigNo*10);
% % plot(ar1totg,'FaceColor','b')
% % plot(ar2totg,'FaceColor','r')
% % plot(ar3totg,'FaceColor','g')

% 
% % areamat(:,length(areamat)+1) = [C01(1); C01(2); totarea];
% 
% % areamat = zeros(length(0:0.5:3),length(0:0.5:6.5));
% 
% areamat(((3e4-(C01(2)))/0.5e4)+1,(C01(1)/0.5e4)+1) = totarea;

function xT = ParametricFormConsensus2(V,J,t,S,x0)
    xT1 = kron(V,eye(2))*(kron(diag(diag(cos(J*t))),eye(2)) + kron(diag(diag(sin(J*t))),S))*kron(inv(V),eye(2));
    xT = xT1*x0;
end



