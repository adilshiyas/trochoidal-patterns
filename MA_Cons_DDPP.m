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

LahR1 = 0.03; LahR2 = 0.03; LahR3 = 0.03;

% dT = 0.0001;
dT = 0.001;
tend = RR; %giving 5 seconds for initial orientation correction
t = 0:dT:tend;
% tcons = 5:dT:RR;

x = zeros(6,length(t));
x2 = x;

u = x;
u2 = x;

SGN = 1;

i = 1;
x(:,i) = x0;
z = kron(D',eye(2))*x(:,i);
u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

x2(:,i) = x0_2;
z2 = kron(D',eye(2))*x2(:,i);
u2(:,i) = -(kron(G*D,eye(2,2)))*z2 + kron((B*D),S)*z2;

for i = 2:length(t)
% for i = 2:length(t)

    x(:,i) = x(:,i-1) + dT*eye(6,6)*(u(:,i-1));
    z = kron(D',eye(2))*x(:,i);
    u(:,i) = -(kron(G*D,eye(2,2)))*z + kron((B*D),S)*z;

    x2(:,i) = x2(:,i-1) + dT*eye(6,6)*(u2(:,i-1));
    z2 = kron(D',eye(2))*x2(:,i);
    u2(:,i) = -(kron(G*D,eye(2,2)))*z2 + kron((B*D),S)*z2;

end

testn = 200;
xR1 = zeros(1,length(t)); yR1 = xR1;
xR2 = xR1; yR2 = xR1; xR3 = xR1; yR3 = xR1;

xR4 = xR1; yR4 = xR1; xR5 = xR1; yR5 = xR1; xR6 = xR1; yR6 = xR1;

diriter = 10;

phInit1 = atan2((x(2,diriter)-x(2,1)),(x(1,diriter)-x(1,1)));
phInit2 = atan2((x(4,diriter)-x(4,1)),(x(3,diriter)-x(3,1)));
phInit3 = atan2((x(6,diriter)-x(6,1)),(x(5,diriter)-x(5,1)));
phInit4 = atan2((x2(2,diriter)-x2(2,1)),(x2(1,diriter)-x2(1,1)));
phInit5 = atan2((x2(4,diriter)-x2(4,1)),(x2(3,diriter)-x2(3,1)));
phInit6 = atan2((x2(6,diriter)-x2(6,1)),(x2(5,diriter)-x2(5,1)));


xInit1 = x(1,1); yInit1 = x(2,1);
xInit2 = x(3,1); yInit2 = x(4,1);
xInit3 = x(5,1); yInit3 = x(6,1);
xInit4 = x2(1,1); yInit4 = x2(2,1);
xInit5 = x2(3,1); yInit5 = x2(4,1);
xInit6 = x2(5,1); yInit6 = x2(6,1);

R1 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R2 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R3 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R4 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R5 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");
R6 = differentialDriveKinematics("TrackWidth", 1, "VehicleInputs","VehicleSpeedHeadingRate");

cont1 = controllerPurePursuit;
cont2 = controllerPurePursuit;
cont3 = controllerPurePursuit;
cont4 = controllerPurePursuit;
cont5 = controllerPurePursuit;
cont6 = controllerPurePursuit;

ppi = 1;

v1 = norm([u(1,ppi) u(2,ppi)]);
v2 = norm([u(3,ppi) u(4,ppi)]);
v3 = norm([u(5,ppi) u(6,ppi)]);
v4 = norm([u2(1,ppi) u2(2,ppi)]);
v5 = norm([u2(3,ppi) u2(4,ppi)]);
v6 = norm([u2(5,ppi) u2(6,ppi)]);

cont1.DesiredLinearVelocity = v1;
cont1.MaxAngularVelocity = 100;
cont1.LookaheadDistance = LahR1;
cont1.Waypoints = [x(1,:)' x(2,:)'];

cont2.DesiredLinearVelocity = v2;
cont2.MaxAngularVelocity = 30;
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

xR1(ppi) = xInit1; yR1(ppi) = yInit1;
xR2(ppi) = xInit2; yR2(ppi) = yInit2;
xR3(ppi) = xInit3; yR3(ppi) = yInit3;
xR4(ppi) = xInit4; yR4(ppi) = yInit4;
xR5(ppi) = xInit5; yR5(ppi) = yInit5;
xR6(ppi) = xInit6; yR6(ppi) = yInit6;

R1CP = [xInit1 yInit1 phInit1]';
R2CP = [xInit2 yInit2 phInit2]';
R3CP = [xInit3 yInit3 phInit3]';
R4CP = [xInit4 yInit4 phInit4]';
R5CP = [xInit5 yInit5 phInit5]';
R6CP = [xInit6 yInit6 phInit6]';

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

   v1 = norm([u(1,ppi) u(2,ppi)]);
   v2 = norm([u(3,ppi) u(4,ppi)]);
   v3 = norm([u(5,ppi) u(6,ppi)]);
   v4 = norm([u2(1,ppi) u2(2,ppi)]);
   v5 = norm([u2(3,ppi) u2(4,ppi)]);
   v6 = norm([u2(5,ppi) u2(6,ppi)]);

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


end



figure(2); clf(2);
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


p1 = plot(xR1(1),yR1(1),'o','MarkerFaceColor','blue');
p2 = plot(xR2(1),yR2(1),'o','MarkerFaceColor','red');
p3 = plot(xR3(1),yR3(1),'o','MarkerFaceColor','green');
p4 = plot(xR4(1),yR4(1),'o','MarkerFaceColor','c');
p5 = plot(xR5(1),yR5(1),'o','MarkerFaceColor','m');
p6 = plot(xR6(1),yR6(1),'o','MarkerFaceColor','k');

for k=2:5:length(xR1)
    p1.XData = xR1(k);
    p1.YData = yR1(k);
    p2.XData = xR2(k);
    p2.YData = yR2(k);
    p3.XData = xR3(k);
    p3.YData = yR3(k);
    p4.XData = xR4(k);
    p4.YData = yR4(k);
    p5.XData = xR5(k);
    p5.YData = yR5(k);
    p6.XData = xR6(k);
    p6.YData = yR6(k);
    drawnow
end


function xT = ParametricFormConsensus2(V,J,t,S,x0)
    xT1 = kron(V,eye(2))*(kron(diag(diag(cos(J*t))),eye(2)) + kron(diag(diag(sin(J*t))),S))*kron(inv(V),eye(2));
    xT = xT1*x0;
end