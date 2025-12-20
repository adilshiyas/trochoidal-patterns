clc;
clear;

k = 2; 
PatternType = 0;
CollDist = 0.5;
CommRange = 10;

%5,12,13 triplet
s1 = 5; s2 = 12; s3 = 13; %5, 12, 13

if(PatternType ==0)
    JR = k + 1;
else
    JR = k - 1;
    JR = -JR;
end

B2 = s2/2;
B3 =  (s3 + s2 - s1 + s1*JR - JR*s2 + JR*s3)/(2*(JR-1));
B1 = B3 - s1;

scaleFactor = 0.1;
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

D = [-1 0;1 -1;0 1];
G = diag([0 0 0]);

L = D*D';
S = [0 -1; 1 0];

[vBL,eBL] = eig(B*L);

C01 = [1600 1000];
C01 = [1610 980];
C01 = [2400 1525];
% C01 = [3000 1800];
C01 = C01*(scaleFactor)^2;

syms x10 x20 y10 y20
PhiA = 0; PhiB = 0;
eqx1 = C01(1) == sqrt((G1*x10 + G2*x20)^2 + (G1*y10 + G2*y20)^2);
eqx2 = C01(2) == sqrt((G1p*x10 + G2p*x20)^2 + (G1p*y10 + G2p*y20)^2);
eqx3 = atan2((G1Phi*y10 + G2Phi*y20),(G1Phi*x10 + G2Phi*x20)) == PhiA;
eqx4 = atan2((G1Phip*y10 + G2Phip*y20),(G1Phip*x10 + G2Phip*x20)) == PhiB;

HardSolve = solve([eqx1,eqx2,eqx3,eqx4],[x10, x20, y10, y20]);
x0 = [round(HardSolve.x10,4);round(HardSolve.y10,4);round(HardSolve.x20,4);round(HardSolve.y20,4);0;0];
x0_Orig = x0;

M = -(kron((G*L),eye(2,2))-kron((B*L),S));
V=null(M');
M1 = [V(:,1)'*(kron(ones(3,1),eye(2,2))); V(:,2)'*(kron(ones(3,1),eye(2,2)))];
M2 = [V(:,1)'*x0; V(:,2)'*x0];
xInf = M1\M2;


x0 = x0 - repmat(xInf,3,1);
x0 = round(x0,3);
M2 = [V(:,1)'*x0; V(:,2)'*x0];
xInf2 = M1\M2;
disp(x0)

eBL = eig(B*L);
reBL = round(diag(eBL),3);
LJ = min(abs(nonzeros(reBL)));
RR = 2*pi/LJ;

x0_2 = ParametricFormConsensus2(vBL,diag(eBL),pi/(2*LJ),S,x0);
x0_2 = round(x0_2,3);
x0_3 = ParametricFormConsensus2(vBL,diag(eBL),pi/(LJ),S,x0);
x0_3 = round(x0_3,3);
x0_4 = ParametricFormConsensus2(vBL,diag(eBL),3*pi/(2*LJ),S,x0);
x0_4 = round(x0_4,3);

LahR = 0.03;
LahR1 = 0.03; LahR2 = 0.03; LahR3 = 0.03;

% dT = 0.0001;
dT = 0.001;
tend = RR; %giving 5 seconds for initial orientation correction
t = 0:dT:tend;
% tcons = 5:dT:RR;

x = zeros(6,length(t));
x2 = x;
x3 = x;
x4 = x;

u = x;
u2 = x;
u3 = x;
u4 = x;

SGN = 1;

i = 1;
x(:,i) = x0;
z = kron(D',eye(2))*x(:,i);
u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

x2(:,i) = x0_2;
z2 = kron(D',eye(2))*x2(:,i);
u2(:,i) = -(kron(G*D,eye(2,2)))*z2 + kron((B*D),S)*z2;

x3(:,i) = x0_3;
z3 = kron(D',eye(2))*x3(:,i);
u3(:,i) = -(kron(G*D,eye(2,2)))*z3 + kron((B*D),S)*z3;

x4(:,i) = x0_4;
z4 = kron(D',eye(2))*x4(:,i);
u4(:,i) = -(kron(G*D,eye(2,2)))*z4 + kron((B*D),S)*z4;

for i = 2:length(t)
% for i = 2:length(t)

    x(:,i) = x(:,i-1) + dT*eye(6,6)*(u(:,i-1));
    z = kron(D',eye(2))*x(:,i);
    u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

    x2(:,i) = x2(:,i-1) + dT*eye(6,6)*(u2(:,i-1));
    z2 = kron(D',eye(2))*x2(:,i);
    u2(:,i) = -(kron(G*D,eye(2,2)))*z2 + kron((B*D),S)*z2;

    x3(:,i) = x3(:,i-1) + dT*eye(6,6)*(u3(:,i-1));
    z3 = kron(D',eye(2))*x3(:,i);
    u3(:,i) = -(kron(G*D,eye(2,2)))*z3 + kron((B*D),S)*z3;

    x4(:,i) = x4(:,i-1) + dT*eye(6,6)*(u4(:,i-1));
    z4 = kron(D',eye(2))*x4(:,i);
    u4(:,i) = -(kron(G*D,eye(2,2)))*z4 + kron((B*D),S)*z4;

end

disp('Cons Complete')

testn = 200;
xR1 = zeros(1,length(t)); yR1 = xR1;
xR2 = xR1; yR2 = xR1; xR3 = xR1; yR3 = xR1;

xR4 = xR1; yR4 = xR1; xR5 = xR1; yR5 = xR1; xR6 = xR1; yR6 = xR1;
xR7 = xR1; yR7 = xR1; xR8 = xR1; yR8 = xR1; xR9 = xR1; yR9 = xR1;
xR10 = xR1; yR10 = xR1; xR11 = xR1; yR11 = xR1; xR12 = xR1; yR12 = xR1;


diriter = 10;

phInit1 = atan2((x(2,diriter)-x(2,1)),(x(1,diriter)-x(1,1)));
phInit2 = atan2((x(4,diriter)-x(4,1)),(x(3,diriter)-x(3,1)));
phInit3 = atan2((x(6,diriter)-x(6,1)),(x(5,diriter)-x(5,1)));
phInit4 = atan2((x2(2,diriter)-x2(2,1)),(x2(1,diriter)-x2(1,1)));
phInit5 = atan2((x2(4,diriter)-x2(4,1)),(x2(3,diriter)-x2(3,1)));
phInit6 = atan2((x2(6,diriter)-x2(6,1)),(x2(5,diriter)-x2(5,1)));
phInit7 = atan2((x3(2,diriter)-x3(2,1)),(x3(1,diriter)-x3(1,1)));
phInit8 = atan2((x3(4,diriter)-x3(4,1)),(x3(3,diriter)-x3(3,1)));
phInit9 = atan2((x3(6,diriter)-x3(6,1)),(x3(5,diriter)-x3(5,1)));
phInit10 = atan2((x4(2,diriter)-x4(2,1)),(x4(1,diriter)-x4(1,1)));
phInit11 = atan2((x4(4,diriter)-x4(4,1)),(x4(3,diriter)-x4(3,1)));
phInit12 = atan2((x4(6,diriter)-x4(6,1)),(x4(5,diriter)-x4(5,1)));


xInit1 = x(1,1); yInit1 = x(2,1);
xInit2 = x(3,1); yInit2 = x(4,1);
xInit3 = x(5,1); yInit3 = x(6,1);
xInit4 = x2(1,1); yInit4 = x2(2,1);
xInit5 = x2(3,1); yInit5 = x2(4,1);
xInit6 = x2(5,1); yInit6 = x2(6,1);
xInit7 = x3(1,1); yInit7 = x3(2,1);
xInit8 = x3(3,1); yInit8 = x3(4,1);
xInit9 = x3(5,1); yInit9 = x3(6,1);
xInit10 = x4(1,1); yInit10 = x4(2,1);
xInit11 = x4(3,1); yInit11 = x4(4,1);
xInit12 = x4(5,1); yInit12 = x4(6,1);

R1 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R2 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R3 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R4 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R5 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R6 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R7 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R8 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R9 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R10 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R11 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R12 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");

cont1 = controllerPurePursuit;
cont2 = controllerPurePursuit;
cont3 = controllerPurePursuit;
cont4 = controllerPurePursuit;
cont5 = controllerPurePursuit;
cont6 = controllerPurePursuit;
cont7 = controllerPurePursuit;
cont8 = controllerPurePursuit;
cont9 = controllerPurePursuit;
cont10 = controllerPurePursuit;
cont11 = controllerPurePursuit;
cont12 = controllerPurePursuit;

ppi = 1;

v1 = norm([u(1,ppi) u(2,ppi)]);
v2 = norm([u(3,ppi) u(4,ppi)]);
v3 = norm([u(5,ppi) u(6,ppi)]);
v4 = norm([u2(1,ppi) u2(2,ppi)]);
v5 = norm([u2(3,ppi) u2(4,ppi)]);
v6 = norm([u2(5,ppi) u2(6,ppi)]);
v7 = norm([u3(1,ppi) u3(2,ppi)]);
v8 = norm([u3(3,ppi) u3(4,ppi)]);
v9 = norm([u3(5,ppi) u3(6,ppi)]);
v10 = norm([u4(1,ppi) u4(2,ppi)]);
v11 = norm([u4(3,ppi) u4(4,ppi)]);
v12 = norm([u4(5,ppi) u4(6,ppi)]);

cont1.DesiredLinearVelocity = v1;
cont1.MaxAngularVelocity = 100;
cont1.LookaheadDistance = LahR1;
cont1.Waypoints = [x(1,:)' x(2,:)'];

cont2.DesiredLinearVelocity = v2;
cont2.MaxAngularVelocity = 100;
cont2.LookaheadDistance = LahR2;
cont2.Waypoints = [x(3,:)' x(4,:)'];

cont3.DesiredLinearVelocity = v3;
cont3.MaxAngularVelocity = 100;
cont3.LookaheadDistance = LahR3;
cont3.Waypoints = [x(5,:)' x(6,:)'];

cont4.DesiredLinearVelocity = v4;
cont4.MaxAngularVelocity = 100;
cont4.LookaheadDistance = LahR1;
cont4.Waypoints = [x2(1,:)' x2(2,:)'];

cont5.DesiredLinearVelocity = v5;
cont5.MaxAngularVelocity = 100;
cont5.LookaheadDistance = LahR2;
cont5.Waypoints = [x2(3,:)' x2(4,:)'];

cont6.DesiredLinearVelocity = v6;
cont6.MaxAngularVelocity = 100;
cont6.LookaheadDistance = LahR3;
cont6.Waypoints = [x2(5,:)' x2(6,:)'];

cont7.DesiredLinearVelocity = v7;
cont7.MaxAngularVelocity = 100;
cont7.LookaheadDistance = LahR;
cont7.Waypoints = [x3(1,:)' x3(2,:)'];

cont8.DesiredLinearVelocity = v8;
cont8.MaxAngularVelocity = 100;
cont8.LookaheadDistance = LahR;
cont8.Waypoints = [x3(3,:)' x3(4,:)'];

cont9.DesiredLinearVelocity = v9;
cont9.MaxAngularVelocity = 100;
cont9.LookaheadDistance = LahR;
cont9.Waypoints = [x3(5,:)' x3(6,:)'];

cont10.DesiredLinearVelocity = v10;
cont10.MaxAngularVelocity = 100;
cont10.LookaheadDistance = LahR;
cont10.Waypoints = [x4(1,:)' x4(2,:)'];

cont11.DesiredLinearVelocity = v11;
cont11.MaxAngularVelocity = 100;
cont11.LookaheadDistance = LahR;
cont11.Waypoints = [x4(3,:)' x4(4,:)'];

cont12.DesiredLinearVelocity = v12;
cont12.MaxAngularVelocity = 100;
cont12.LookaheadDistance = LahR;
cont12.Waypoints = [x4(5,:)' x4(6,:)'];

xR1(ppi) = xInit1; yR1(ppi) = yInit1;
xR2(ppi) = xInit2; yR2(ppi) = yInit2;
xR3(ppi) = xInit3; yR3(ppi) = yInit3;
xR4(ppi) = xInit4; yR4(ppi) = yInit4;
xR5(ppi) = xInit5; yR5(ppi) = yInit5;
xR6(ppi) = xInit6; yR6(ppi) = yInit6;
xR7(ppi) = xInit7; yR7(ppi) = yInit7;
xR8(ppi) = xInit8; yR8(ppi) = yInit8;
xR9(ppi) = xInit9; yR9(ppi) = yInit9;
xR10(ppi) = xInit10; yR10(ppi) = yInit10;
xR11(ppi) = xInit11; yR11(ppi) = yInit11;
xR12(ppi) = xInit12; yR12(ppi) = yInit12;

R1CP = [xInit1 yInit1 phInit1]';
R2CP = [xInit2 yInit2 phInit2]';
R3CP = [xInit3 yInit3 phInit3]';
R4CP = [xInit4 yInit4 phInit4]';
R5CP = [xInit5 yInit5 phInit5]';
R6CP = [xInit6 yInit6 phInit6]';
R7CP = [xInit7 yInit7 phInit7]';
R8CP = [xInit8 yInit8 phInit8]';
R9CP = [xInit9 yInit9 phInit9]';
R10CP = [xInit10 yInit10 phInit10]';
R11CP = [xInit11 yInit11 phInit11]';
R12CP = [xInit12 yInit12 phInit12]';


[v_1,omega_1] = cont1(R1CP);
vel1 = derivative(R1,R1CP,[v_1,omega_1]);
[v_2,omega_2] = cont2(R2CP);
vel2 = derivative(R2,R2CP,[v_2,omega_2]);
[v_3,omega_3] = cont3(R3CP);
vel3 = derivative(R3,R3CP,[v_3,omega_3]);
[v_4,omega_4] = cont4(R4CP);
vel4 = derivative(R4,R4CP,[v_4,omega_4]);
[v_5,omega_5] = cont5(R5CP);
vel5 = derivative(R5,R5CP,[v_5,omega_5]);
[v_6,omega_6] = cont6(R6CP);
vel6 = derivative(R6,R6CP,[v_6,omega_6]);
[v_7,omega_7] = cont7(R7CP);
vel7 = derivative(R7,R7CP,[v_7,omega_7]);
[v_8,omega_8] = cont8(R8CP);
vel8 = derivative(R8,R8CP,[v_8,omega_8]);
[v_9,omega_9] = cont9(R9CP);
vel9 = derivative(R9,R9CP,[v_9,omega_9]);
[v_10,omega_10] = cont10(R10CP);
vel10 = derivative(R10,R10CP,[v_10,omega_10]);
[v_11,omega_11] = cont11(R11CP);
vel11 = derivative(R11,R11CP,[v_11,omega_11]);
[v_12,omega_12] = cont12(R12CP);
vel12 = derivative(R12,R12CP,[v_12,omega_12]);

for ppi = 2:length(t)

   R1CP = R1CP + vel1*dT;
   xR1(ppi) = R1CP(1); yR1(ppi) = R1CP(2);
   R2CP = R2CP + vel2*dT;
   xR2(ppi) = R2CP(1); yR2(ppi) = R2CP(2);
   R3CP = R3CP + vel3*dT;
   xR3(ppi) = R3CP(1); yR3(ppi) = R3CP(2);
   R4CP = R4CP + vel4*dT;
   xR4(ppi) = R4CP(1); yR4(ppi) = R4CP(2);
   R5CP = R5CP + vel5*dT;
   xR5(ppi) = R5CP(1); yR5(ppi) = R5CP(2);
   R6CP = R6CP + vel6*dT;
   xR6(ppi) = R6CP(1); yR6(ppi) = R6CP(2);
   R7CP = R7CP + vel7*dT;
   xR7(ppi) = R7CP(1); yR7(ppi) = R7CP(2);
   R8CP = R8CP + vel8*dT;
   xR8(ppi) = R8CP(1); yR8(ppi) = R8CP(2);
   R9CP = R9CP + vel9*dT;
   xR9(ppi) = R9CP(1); yR9(ppi) = R9CP(2);
   R10CP = R10CP + vel10*dT;
   xR10(ppi) = R10CP(1); yR10(ppi) = R10CP(2);
   R11CP = R11CP + vel11*dT;
   xR11(ppi) = R11CP(1); yR11(ppi) = R11CP(2);
   R12CP = R12CP + vel12*dT;
   xR12(ppi) = R12CP(1); yR12(ppi) = R12CP(2);

   v1 = norm([u(1,ppi) u(2,ppi)]);
   v2 = norm([u(3,ppi) u(4,ppi)]);
   v3 = norm([u(5,ppi) u(6,ppi)]);
   v4 = norm([u2(1,ppi) u2(2,ppi)]);
   v5 = norm([u2(3,ppi) u2(4,ppi)]);
   v6 = norm([u2(5,ppi) u2(6,ppi)]);
   v7 = norm([u3(1,ppi) u3(2,ppi)]);
   v8 = norm([u3(3,ppi) u3(4,ppi)]);
   v9 = norm([u3(5,ppi) u3(6,ppi)]);
   v10 = norm([u4(1,ppi) u4(2,ppi)]);
   v11 = norm([u4(3,ppi) u4(4,ppi)]);
   v12 = norm([u4(5,ppi) u4(6,ppi)]);

   cont1.DesiredLinearVelocity = v1;
   [v_1,omega_1] = cont1(R1CP);
   vel1 = derivative(R1,R1CP,[v_1,omega_1]);

   cont2.DesiredLinearVelocity = v2;
   [v_2,omega_2] = cont2(R2CP);
   vel2 = derivative(R2,R2CP,[v_2,omega_2]);

   cont3.DesiredLinearVelocity = v3;
   [v_3,omega_3] = cont3(R3CP);
   vel3 = derivative(R3,R3CP,[v_3,omega_3]);
  
   cont4.DesiredLinearVelocity = v4;
   [v_4,omega_4] = cont4(R4CP);
   vel4 = derivative(R4,R4CP,[v_4,omega_4]);

   cont5.DesiredLinearVelocity = v5;
   [v_5,omega_5] = cont5(R5CP);
   vel5 = derivative(R5,R5CP,[v_5,omega_5]);

   cont6.DesiredLinearVelocity = v6;
   [v_6,omega_6] = cont6(R6CP);
   vel6 = derivative(R6,R6CP,[v_6,omega_6]);

   cont7.DesiredLinearVelocity = v7;
   [v_7,omega_7] = cont7(R7CP);
   vel7 = derivative(R7,R7CP,[v_7,omega_7]);

   cont8.DesiredLinearVelocity = v8;
   [v_8,omega_8] = cont8(R8CP);
   vel8 = derivative(R8,R8CP,[v_8,omega_8]);

   cont9.DesiredLinearVelocity = v9;
   [v_9,omega_9] = cont9(R9CP);
   vel9 = derivative(R9,R9CP,[v_9,omega_9]);

   cont10.DesiredLinearVelocity = v10;
   [v_10,omega_10] = cont10(R10CP);
   vel10 = derivative(R10,R10CP,[v_10,omega_10]);

   cont11.DesiredLinearVelocity = v11;
   [v_11,omega_11] = cont11(R11CP);
   vel11 = derivative(R11,R11CP,[v_11,omega_11]);

   cont12.DesiredLinearVelocity = v12;
   [v_12,omega_12] = cont12(R12CP);
   vel12 = derivative(R12,R12CP,[v_12,omega_12]);

end



% figure(2); clf(2);
% plot(x(1,:),x(2,:),'b')
% hold on
% plot(x(3,:),x(4,:),'r')
% plot(x(5,:),x(6,:),'g')
% plot(x(1,1),x(2,1),'ob')
% plot(x(3,1),x(4,1),'or')
% plot(x(5,1),x(6,1),'og')
% 
% % tht  = 0:0.1:2*pi;
% % plot(x(1,1)+LahR1*cos(tht),x(2,1)+LahR1*sin(tht),'b')
% % plot(x(3,1)+LahR2*cos(tht),x(4,1)+LahR2*sin(tht),'r')
% % plot(x(5,1)+LahR3*cos(tht),x(6,1)+LahR3*sin(tht),'g')
% 
% 
% p1 = plot(xR1(1),yR1(1),'o','MarkerFaceColor','blue');
% p2 = plot(xR2(1),yR2(1),'o','MarkerFaceColor','red');
% p3 = plot(xR3(1),yR3(1),'o','MarkerFaceColor','green');
% p4 = plot(xR4(1),yR4(1),'o','MarkerFaceColor','c');
% p5 = plot(xR5(1),yR5(1),'o','MarkerFaceColor','m');
% p6 = plot(xR6(1),yR6(1),'o','MarkerFaceColor','k');
% p7 = plot(xR7(1),yR7(1),'o','MarkerFaceColor','[0 0.4470 0.7410]');
% p8 = plot(xR8(1),yR8(1),'o','MarkerFaceColor','[0.8500 0.3250 0.0980]');
% p9 = plot(xR9(1),yR9(1),'o','MarkerFaceColor','[0.4660 0.6740 0.1880]');
% p10 = plot(xR10(1),yR10(1),'o','MarkerFaceColor','[0.3010 0.7450 0.9330]');
% p11 = plot(xR11(1),yR11(1),'o','MarkerFaceColor','[0.6350 0.0780 0.1840]');
% p12 = plot(xR12(1),yR12(1),'o','MarkerFaceColor','#013220');
% 
% for k=2:5:length(xR1)
%     p1.XData = xR1(k);
%     p1.YData = yR1(k);
%     p2.XData = xR2(k);
%     p2.YData = yR2(k);
%     p3.XData = xR3(k);
%     p3.YData = yR3(k);
%     p4.XData = xR4(k);
%     p4.YData = yR4(k);
%     p5.XData = xR5(k);
%     p5.YData = yR5(k);
%     p6.XData = xR6(k);
%     p6.YData = yR6(k);
%     p7.XData = xR7(k);
%     p7.YData = yR7(k);
%     p8.XData = xR8(k);
%     p8.YData = yR8(k);
%     p9.XData = xR9(k);
%     p9.YData = yR9(k);
%     p10.XData = xR10(k);
%     p10.YData = yR10(k);
%     p11.XData = xR11(k);
%     p11.YData = yR11(k);
%     p12.XData = xR12(k);
%     p12.YData = yR12(k);
%     drawnow
% end

figure(3); clf(3);
% figure('Position',[0 0 1500 1000]);
plot(x(1,:),x(2,:),'b')
hold on
plot(x(3,:),x(4,:),'r')
plot(x(5,:),x(6,:),'g')
plot(x(1,1),x(2,1),'ob')
plot(x(3,1),x(4,1),'or')
plot(x(5,1),x(6,1),'og')

% tht  = 0:0.1:2*pi;
% plot(x(1,1)+LahR1*cos(tht),x(2,1)+LahR1*sin(tht),'b')
% plot(x(3,1)+LahR2*cos(tht),x(4,1)+LahR2*sin(tht),'r')
% plot(x(5,1)+LahR3*cos(tht),x(6,1)+LahR3*sin(tht),'g')

thc = 0:0.1:2*pi;
rbot = 0.45/2;

op1 = plot(polyshape(xR1(1)+ rbot*cos(thc),yR1(1)+ rbot*sin(thc))); op1.FaceColor = 'b';
op2 = plot(polyshape(xR2(1)+ rbot*cos(thc),yR2(1)+ rbot*sin(thc))); op2.FaceColor = 'r';
op3 = plot(polyshape(xR3(1)+ rbot*cos(thc),yR3(1)+ rbot*sin(thc))); op3.FaceColor = 'g';
op4 = plot(polyshape(xR4(1)+ rbot*cos(thc),yR4(1)+ rbot*sin(thc))); op4.FaceColor = 'c';
op5 = plot(polyshape(xR5(1)+ rbot*cos(thc),yR5(1)+ rbot*sin(thc))); op5.FaceColor = 'm';
op6 = plot(polyshape(xR6(1)+ rbot*cos(thc),yR6(1)+ rbot*sin(thc))); op6.FaceColor = 'k';
op7 = plot(polyshape(xR7(1)+ rbot*cos(thc),yR7(1)+ rbot*sin(thc))); op7.FaceColor = '[0 0.4470 0.7410]';
op8 = plot(polyshape(xR8(1)+ rbot*cos(thc),yR8(1)+ rbot*sin(thc))); op8.FaceColor = '[0.8500 0.3250 0.0980]';
op9 = plot(polyshape(xR9(1)+ rbot*cos(thc),yR9(1)+ rbot*sin(thc))); op9.FaceColor = '[0.4660 0.6740 0.1880]';
op10 = plot(polyshape(xR10(1)+ rbot*cos(thc),yR10(1)+ rbot*sin(thc))); op10.FaceColor = '[0.3010 0.7450 0.9330]';
op11 = plot(polyshape(xR11(1)+ rbot*cos(thc),yR11(1)+ rbot*sin(thc))); op11.FaceColor = '[0.6350 0.0780 0.1840]';
op12 = plot(polyshape(xR12(1)+ rbot*cos(thc),yR12(1)+ rbot*sin(thc))); op12.FaceColor = '#013220';


p1 = plot(xR1(1),yR1(1),'o','MarkerFaceColor','blue');
p2 = plot(xR2(1),yR2(1),'o','MarkerFaceColor','red');
p3 = plot(xR3(1),yR3(1),'o','MarkerFaceColor','green');
p4 = plot(xR4(1),yR4(1),'o','MarkerFaceColor','c');
p5 = plot(xR5(1),yR5(1),'o','MarkerFaceColor','m');
p6 = plot(xR6(1),yR6(1),'o','MarkerFaceColor','k');
p7 = plot(xR7(1),yR7(1),'o','MarkerFaceColor','[0 0.4470 0.7410]');
p8 = plot(xR8(1),yR8(1),'o','MarkerFaceColor','[0.8500 0.3250 0.0980]');
p9 = plot(xR9(1),yR9(1),'o','MarkerFaceColor','[0.4660 0.6740 0.1880]');
p10 = plot(xR10(1),yR10(1),'o','MarkerFaceColor','[0.3010 0.7450 0.9330]');
p11 = plot(xR11(1),yR11(1),'o','MarkerFaceColor','[0.6350 0.0780 0.1840]');
p12 = plot(xR12(1),yR12(1),'o','MarkerFaceColor','#013220');
ax = gca;
ax.PlotBoxAspectRatio = [3 2 1];

% for k=2:10:length(xR1)
%     delete(op1)
%     delete(op2)
%     delete(op3)
%     delete(op4)
%     delete(op5)
%     delete(op6)
%     delete(op7)
%     delete(op8)
%     delete(op9)
%     delete(op10)
%     delete(op11)
%     delete(op12)
%     
% 
% 
%     p1.XData = xR1(k);
%     p1.YData = yR1(k);
%     p2.XData = xR2(k);
%     p2.YData = yR2(k);
%     p3.XData = xR3(k);
%     p3.YData = yR3(k);
%     p4.XData = xR4(k);
%     p4.YData = yR4(k);
%     p5.XData = xR5(k);
%     p5.YData = yR5(k);
%     p6.XData = xR6(k);
%     p6.YData = yR6(k);
%     p7.XData = xR7(k);
%     p7.YData = yR7(k);
%     p8.XData = xR8(k);
%     p8.YData = yR8(k);
%     p9.XData = xR9(k);
%     p9.YData = yR9(k);
%     p10.XData = xR10(k);
%     p10.YData = yR10(k);
%     p11.XData = xR11(k);
%     p11.YData = yR11(k);
%     p12.XData = xR12(k);
%     p12.YData = yR12(k);
% 
%     op1 = plot(polyshape(xR1(k)+ rbot*cos(thc),yR1(k)+ rbot*sin(thc))); op1.FaceColor = 'b';
%     op2 = plot(polyshape(xR2(k)+ rbot*cos(thc),yR2(k)+ rbot*sin(thc))); op2.FaceColor = 'r';
%     op3 = plot(polyshape(xR3(k)+ rbot*cos(thc),yR3(k)+ rbot*sin(thc))); op3.FaceColor = 'g';
%     op4 = plot(polyshape(xR4(k)+ rbot*cos(thc),yR4(k)+ rbot*sin(thc))); op4.FaceColor = 'c';
%     op5 = plot(polyshape(xR5(k)+ rbot*cos(thc),yR5(k)+ rbot*sin(thc))); op5.FaceColor = 'm';
%     op6 = plot(polyshape(xR6(k)+ rbot*cos(thc),yR6(k)+ rbot*sin(thc))); op6.FaceColor = 'k';
%     op7 = plot(polyshape(xR7(k)+ rbot*cos(thc),yR7(k)+ rbot*sin(thc))); op7.FaceColor = '[0 0.4470 0.7410]';
%     op8 = plot(polyshape(xR8(k)+ rbot*cos(thc),yR8(k)+ rbot*sin(thc))); op8.FaceColor = '[0.8500 0.3250 0.0980]';
%     op9 = plot(polyshape(xR9(k)+ rbot*cos(thc),yR9(k)+ rbot*sin(thc))); op9.FaceColor = '[0.4660 0.6740 0.1880]';
%     op10 = plot(polyshape(xR10(k)+ rbot*cos(thc),yR10(k)+ rbot*sin(thc))); op10.FaceColor = '[0.3010 0.7450 0.9330]';
%     op11 = plot(polyshape(xR11(k)+ rbot*cos(thc),yR11(k)+ rbot*sin(thc))); op11.FaceColor = '[0.6350 0.0780 0.1840]';
%     op12 = plot(polyshape(xR12(k)+ rbot*cos(thc),yR12(k)+ rbot*sin(thc))); op12.FaceColor = '#013220';
%     drawnow
% end

function xT = ParametricFormConsensus2(V,J,t,S,x0)
    xT1 = kron(V,eye(2))*(kron(diag(diag(cos(J*t))),eye(2)) + kron(diag(diag(sin(J*t))),S))*kron(inv(V),eye(2));
    xT = xT1*x0;
end