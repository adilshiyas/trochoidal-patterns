

k = 2;
PatternType = 0; %0 = Epi, 1 = Hypo
D = [-1 0;1 -1;0 1];
G = diag([0 0 0]);

L = D*D';
S = [0 -1; 1 0];

%Triple choices
% s1 = 3; s2 = 4; s3 = 5;           %3, 4, 5
s1 = 5; s2 = 12; s3 = 13;         %5, 12, 13
% s1 = s1/10; s2 = s2/10; s3 = s3/10;
% s1 = 8; s2 = 15; s3 = 17;         %8, 15, 17
% s1 = 7; s2 = 24; s3 = 25;         %7, 24, 25
% s1 = 28; s2 = 45; s3 = 53;         %28, 45, 53
% s1 = 20; s2 = 21; s3 = 29;         %20, 21, 29
% s1 = 12; s2 = 35; s3 = 37;         %12, 35, 37
% s1 = 9; s2 = 40; s3 = 41;         %9, 40, 41
% s1 = 11; s2 = 60; s3 = 61;         %11, 60, 61
% s1 = 16; s2 = 63; s3 = 65;         %16, 63, 65

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


% scaleFactor = 1/(10*1.9365); %polyshape down to 0-1 by 0-1
scaleFactor = 1;
B1 = B1*scaleFactor;
B2 = B2*scaleFactor;
B3 = B3*scaleFactor;

B = diag([B1 B2 B3]);
[vBL, eBL] = eig(B*L);

eBL_zind = find(diag(eBL)==0);
eBL_maxind = find(diag(eBL)==max(diag(eBL)));
eBL_minind = find(diag(eBL)==min(diag(eBL)));
invvBL = inv(vBL);

a = B1/2 + B2 + B3/2;
b = (B1^2 - 2*B1*B3 + 4*B2^2 + B3^2)^(1/2)/2;
alpha3r = B3/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha3d = alpha3r;
alpha2r = (B3 - (a-b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha2d = (B3 - (a+b))/(2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1r = (B1*(a-b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));
alpha1d = (B1*(a+b) - B1*B2 - B1*B3)/(2*B2*(B1*B2*2*b + B2*B3*2*b + B1*B3*2*b));

alphamat = [alpha1r alpha1d; alpha2r alpha2d; alpha3r alpha3d];

G2 = 2*b*B1 + B1^2 - B1*B3 - 2*B2*B3;
G3 = -(2*b*B1 + B1^2 - B1*B3 -B2*B3 + 2*b*B2 + 2*B2^2 + B1*B2);
G1 = B2*B3 + 2*b*B2 + 2*B2^2 +B1*B2;

G3p = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2p = -(B1*2*b - B1^2 + B1*B3 + 2*B2*B3);
G1p = B1*B2 + 2*B2^2 + B2*B3 - B2*2*b;


G3Phi = B1*2*b + B2*2*b + B1^2 + 2*B2^2 + B1*B2 - B1*B3 - B2*B3;
G2Phi = -(B1*2*b + B1^2 - B1*B3 - 2*B2*B3);
G1Phi = -(B2*B1 + 2*B2^2 + B2*B3 + B2*2*b);

G3Phip = B1*2*b + B2*2*b - B1^2 - 2*B2^2 - B1*B2 + B1*B3 + B2*B3;
G2Phip = -(B1*2*b -B1^2 + B1*B3 + 2*B2*B3);
G1Phip = B2*B1 + 2*B2^2 + B3*B2 - B2*2*b;

% Constraint values
CValLim = 100000;
dmindelt = 0;
Dmin = 0.2;
Dmax = 2;
CollThresh = 0.35;
RangeThresh = 4;

% abs(alpha1r)*Rc - abs(alpha1d)*dc >= Dmin
CValLim = CValLim*scaleFactor;
ttestno = 3;

iPhiA = 0;
iPhiB = 0;

% Finding solutions to Dmin linear inequalities
syms dc Rc
% eq = abs(alpha1r)*Rc - abs(alpha1d)*dc == Dmin;
% eqneg = abs(alpha1r)*Rc - abs(alpha1d)*dc == -Dmin;
pts1 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha1r)*Rc - abs(alpha1d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha1r)*Rc - abs(alpha1d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV;      
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha1r)*Rc - abs(alpha1d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha1r)*Rc - abs(alpha1d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha2r)*Rc - abs(alpha2d)*dc == Dmin;
% eqneg = abs(alpha2r)*Rc - abs(alpha2d)*dc == -Dmin;
pts2 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha2r)*Rc - abs(alpha2d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha2r)*Rc - abs(alpha2d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha2r)*Rc - abs(alpha2d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha2r)*Rc - abs(alpha2d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha3r)*Rc - abs(alpha3d)*dc == Dmin;
% eqneg = abs(alpha3r)*Rc - abs(alpha3d)*dc == -Dmin;
pts3 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha3r)*Rc - abs(alpha3d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha3r)*Rc - abs(alpha3d)*dc == Dmin;
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha3r)*Rc - abs(alpha3d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha3r)*Rc - abs(alpha3d)*dc == -Dmin;
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV;  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1(1,:) == min(pts1(1,:)));
w = find(pts1(2,q) == max(pts1(2,q)));
e = find(pts1(2,q) == min(pts1(2,q)));
r = find(pts1(1,:) == max(pts1(1,:)));
t = find(pts1(2,r) == min(pts1(2,r)));
y = find(pts1(2,r) == max(pts1(2,r)));

pts1t = pts1(:,q(w));
pts1t(:,end+1) = pts1(:,q(e));
pts1t(:,end+1) = pts1(:,r(t));
pts1t(:,end+1) = pts1(:,r(y));
pts1 = pts1t;

q = find(pts2(1,:) == min(pts2(1,:)));
w = find(pts2(2,q) == max(pts2(2,q)));
e = find(pts2(2,q) == min(pts2(2,q)));
r = find(pts2(1,:) == max(pts2(1,:)));
t = find(pts2(2,r) == min(pts2(2,r)));
y = find(pts2(2,r) == max(pts2(2,r)));

pts2t = pts2(:,q(w));
pts2t(:,end+1) = pts2(:,q(e));
pts2t(:,end+1) = pts2(:,r(t));
pts2t(:,end+1) = pts2(:,r(y));
pts2 = pts2t;

q = find(pts3(1,:) == min(pts3(1,:)));
w = find(pts3(2,q) == max(pts3(2,q)));
e = find(pts3(2,q) == min(pts3(2,q)));
r = find(pts3(1,:) == max(pts3(1,:)));
t = find(pts3(2,r) == min(pts3(2,r)));
y = find(pts3(2,r) == max(pts3(2,r)));

pts3t = pts3(:,q(w));
pts3t(:,end+1) = pts3(:,q(e));
pts3t(:,end+1) = pts3(:,r(t));
pts3t(:,end+1) = pts3(:,r(y));
pts3 = pts3t;


pgonNoguard1 = polyshape([pts1(1,1),pts1(1,2),pts1(1,3),pts1(1,4)],[pts1(2,1),pts1(2,2),pts1(2,3),pts1(2,4)]);
pgonNoguard2 = polyshape([pts2(1,1),pts2(1,2),pts2(1,3),pts2(1,4)],[pts2(2,1),pts2(2,2),pts2(2,3),pts2(2,4)]);
pgonNoguard3 = polyshape([pts3(1,1),pts3(1,2),pts3(1,3),pts3(1,4)],[pts3(2,1),pts3(2,2),pts3(2,3),pts3(2,4)]);
pgonNoguard = union(pgonNoguard1,pgonNoguard2);
pgonNoguard = union(pgonNoguard,pgonNoguard3);
pgonGuard = subtract(pgon,pgonNoguard);


% For ensuring R+r+d doesnt cross a certain value (Dmax)
% Finding solutions to Dmax linear inequalities

syms dc Rc
% eq = abs(alpha1r)*Rc + abs(alpha1d)*dc == Dmax;
% eqneg = abs(alpha1r)*Rc + abs(alpha1d)*dc == -Dmax;
pts1m = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha1r)*Rc + abs(alpha1d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts1m(1,i) = Rc;
            pts1m(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha1r)*Rc + abs(alpha1d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts1m(1,i) = Rc;
            pts1m(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha1r)*Rc + abs(alpha1d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts1m(1,i) = Rc;
            pts1m(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha1r)*Rc + abs(alpha1d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts1m(1,i) = Rc;
            pts1m(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha2r)*Rc + abs(alpha2d)*dc == Dmax;
% eqneg = abs(alpha2r)*Rc + abs(alpha2d)*dc == -Dmax;
pts2m = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha2r)*Rc + abs(alpha2d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts2m(1,i) = Rc;
            pts2m(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha2r)*Rc + abs(alpha2d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts2m(1,i) = Rc;
            pts2m(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha2r)*Rc + abs(alpha2d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts2m(1,i) = Rc;
            pts2m(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha2r)*Rc + abs(alpha2d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts2m(1,i) = Rc;
            pts2m(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha3r)*Rc + abs(alpha3d)*dc == Dmax;
% eqneg = abs(alpha3r)*Rc + abs(alpha3d)*dc == -Dmax;
pts3m = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha3r)*Rc + abs(alpha3d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts3m(1,i) = Rc;
            pts3m(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha3r)*Rc + abs(alpha3d)*dc == Dmax;
            SLV = solve(eq,dc);
            pts3m(1,i) = Rc;
            pts3m(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha3r)*Rc + abs(alpha3d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts3m(1,i) = Rc;
            pts3m(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha3r)*Rc + abs(alpha3d)*dc == -Dmax;
            SLV = solve(eqneg,dc);
            pts3m(1,i) = Rc;
            pts3m(2,i) = SLV;  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1m(1,:) == min(pts1m(1,:)));
w = find(pts1m(2,q) == max(pts1m(2,q)));
e = find(pts1m(2,q) == min(pts1m(2,q)));
r = find(pts1m(1,:) == max(pts1m(1,:)));
t = find(pts1m(2,r) == min(pts1m(2,r)));
y = find(pts1m(2,r) == max(pts1m(2,r)));

pts1mt = pts1m(:,q(w));
pts1mt(:,end+1) = pts1m(:,q(e));
pts1mt(:,end+1) = pts1m(:,r(t));
pts1mt(:,end+1) = pts1m(:,r(y));
pts1m = pts1mt;

q = find(pts2m(1,:) == min(pts2m(1,:)));
w = find(pts2m(2,q) == max(pts2m(2,q)));
e = find(pts2m(2,q) == min(pts2m(2,q)));
r = find(pts2m(1,:) == max(pts2m(1,:)));
t = find(pts2m(2,r) == min(pts2m(2,r)));
y = find(pts2m(2,r) == max(pts2m(2,r)));

pts2mt = pts2m(:,q(w));
pts2mt(:,end+1) = pts2m(:,q(e));
pts2mt(:,end+1) = pts2m(:,r(t));
pts2mt(:,end+1) = pts2m(:,r(y));
pts2m = pts2mt;

q = find(pts3m(1,:) == min(pts3m(1,:)));
w = find(pts3m(2,q) == max(pts3m(2,q)));
e = find(pts3m(2,q) == min(pts3m(2,q)));
r = find(pts3m(1,:) == max(pts3m(1,:)));
t = find(pts3m(2,r) == min(pts3m(2,r)));
y = find(pts3m(2,r) == max(pts3m(2,r)));

pts3mt = pts3m(:,q(w));
pts3mt(:,end+1) = pts3m(:,q(e));
pts3mt(:,end+1) = pts3m(:,r(t));
pts3mt(:,end+1) = pts3m(:,r(y));
pts3m = pts3mt;


pgonMaxDist1 = polyshape([pts1m(1,1),pts1m(1,2),pts1m(1,3),pts1m(1,4)],[pts1m(2,1),pts1m(2,2),pts1m(2,3),pts1m(2,4)]);
pgonMaxDist2 = polyshape([pts2m(1,1),pts2m(1,2),pts2m(1,3),pts2m(1,4)],[pts2m(2,1),pts2m(2,2),pts2m(2,3),pts2m(2,4)]);
pgonMaxDist3 = polyshape([pts3m(1,1),pts3m(1,2),pts3m(1,3),pts3m(1,4)],[pts3m(2,1),pts3m(2,2),pts3m(2,3),pts3m(2,4)]);
pgonMaxDist = intersect(pgonMaxDist1,pgonMaxDist2);
pgonMaxDist = intersect(pgonMaxDist,pgonMaxDist3);
% pgonGuard = subtract(pgon,pgonNoguard);

% Ensure that any 2 agents do not collide (come within CollThresh distance
% of each other
% Finding solutions to DCollision linear inequalities

syms dc Rc
% eq = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == CollThresh;
% eqneg = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == -CollThresh;
pts1c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha1r-alpha2r)*Rc - abs(alpha1d-alpha2d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;  
        end

    end
end

eq = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == CollThresh;
eqneg = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == -CollThresh;
pts2c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha2r-alpha3r)*Rc - abs(alpha2d-alpha3d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == CollThresh;
% eqneg = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == -CollThresh;
pts3c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == CollThresh;
            SLV = solve(eq,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha3r-alpha1r)*Rc - abs(alpha3d-alpha1d)*dc == -CollThresh;
            SLV = solve(eqneg,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1c(1,:) == min(pts1c(1,:)));
w = find(pts1c(2,q) == max(pts1c(2,q)));
e = find(pts1c(2,q) == min(pts1c(2,q)));
r = find(pts1c(1,:) == max(pts1c(1,:)));
t = find(pts1c(2,r) == min(pts1c(2,r)));
y = find(pts1c(2,r) == max(pts1c(2,r)));

pts1ct = pts1c(:,q(w));
pts1ct(:,end+1) = pts1c(:,q(e));
pts1ct(:,end+1) = pts1c(:,r(t));
pts1ct(:,end+1) = pts1c(:,r(y));
pts1c = pts1ct;

q = find(pts2c(1,:) == min(pts2c(1,:)));
w = find(pts2c(2,q) == max(pts2c(2,q)));
e = find(pts2c(2,q) == min(pts2c(2,q)));
r = find(pts2c(1,:) == max(pts2c(1,:)));
t = find(pts2c(2,r) == min(pts2c(2,r)));
y = find(pts2c(2,r) == max(pts2c(2,r)));

pts2ct = pts2c(:,q(w));
pts2ct(:,end+1) = pts2c(:,q(e));
pts2ct(:,end+1) = pts2c(:,r(t));
pts2ct(:,end+1) = pts2c(:,r(y));
pts2c = pts2ct;

q = find(pts3c(1,:) == min(pts3c(1,:)));
w = find(pts3c(2,q) == max(pts3c(2,q)));
e = find(pts3c(2,q) == min(pts3c(2,q)));
r = find(pts3c(1,:) == max(pts3c(1,:)));
t = find(pts3c(2,r) == min(pts3c(2,r)));
y = find(pts3c(2,r) == max(pts3c(2,r)));

pts3ct = pts3c(:,q(w));
pts3ct(:,end+1) = pts3c(:,q(e));
pts3ct(:,end+1) = pts3c(:,r(t));
pts3ct(:,end+1) = pts3c(:,r(y));
pts3c = pts3ct;


pgonColl12 = polyshape([pts1c(1,1),pts1c(1,2),pts1c(1,3),pts1c(1,4)],[pts1c(2,1),pts1c(2,2),pts1c(2,3),pts1c(2,4)]);
pgonColl23 = polyshape([pts2c(1,1),pts2c(1,2),pts2c(1,3),pts2c(1,4)],[pts2c(2,1),pts2c(2,2),pts2c(2,3),pts2c(2,4)]);
pgonColl31 = polyshape([pts3c(1,1),pts3c(1,2),pts3c(1,3),pts3c(1,4)],[pts3c(2,1),pts3c(2,2),pts3c(2,3),pts3c(2,4)]);
pgonColl = union(pgonColl12,pgonColl23);
pgonColl = union(pgonColl,pgonColl31);
pgonNoColl = subtract(pgon,pgonColl);

% Ensure that any 2 agents do not collide (come within CollThresh distance
% of each other
% Finding solutions to DRange linear inequalities

syms dc Rc
% eq = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == RangeThresh;
% eqneg = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == -RangeThresh;
pts1c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha1r-alpha2r)*Rc + abs(alpha1d-alpha2d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts1c(1,i) = Rc;
            pts1c(2,i) = SLV;  
        end

    end
end

eq = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == RangeThresh;
eqneg = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == -RangeThresh;
pts2c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha2r-alpha3r)*Rc + abs(alpha2d-alpha3d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts2c(1,i) = Rc;
            pts2c(2,i) = SLV;  
        end

    end
end

% eq = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == RangeThresh;
% eqneg = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == -RangeThresh;
pts3c = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;
        else
            Rc = -CValLim;
            eq = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == RangeThresh;
            SLV = solve(eq,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;
        else
            Rc = -CValLim;
            eqneg = abs(alpha3r-alpha1r)*Rc + abs(alpha3d-alpha1d)*dc == -RangeThresh;
            SLV = solve(eqneg,dc);
            pts3c(1,i) = Rc;
            pts3c(2,i) = SLV;  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1c(1,:) == min(pts1c(1,:)));
w = find(pts1c(2,q) == max(pts1c(2,q)));
e = find(pts1c(2,q) == min(pts1c(2,q)));
r = find(pts1c(1,:) == max(pts1c(1,:)));
t = find(pts1c(2,r) == min(pts1c(2,r)));
y = find(pts1c(2,r) == max(pts1c(2,r)));

pts1ct = pts1c(:,q(w));
pts1ct(:,end+1) = pts1c(:,q(e));
pts1ct(:,end+1) = pts1c(:,r(t));
pts1ct(:,end+1) = pts1c(:,r(y));
pts1c = pts1ct;

q = find(pts2c(1,:) == min(pts2c(1,:)));
w = find(pts2c(2,q) == max(pts2c(2,q)));
e = find(pts2c(2,q) == min(pts2c(2,q)));
r = find(pts2c(1,:) == max(pts2c(1,:)));
t = find(pts2c(2,r) == min(pts2c(2,r)));
y = find(pts2c(2,r) == max(pts2c(2,r)));

pts2ct = pts2c(:,q(w));
pts2ct(:,end+1) = pts2c(:,q(e));
pts2ct(:,end+1) = pts2c(:,r(t));
pts2ct(:,end+1) = pts2c(:,r(y));
pts2c = pts2ct;

q = find(pts3c(1,:) == min(pts3c(1,:)));
w = find(pts3c(2,q) == max(pts3c(2,q)));
e = find(pts3c(2,q) == min(pts3c(2,q)));
r = find(pts3c(1,:) == max(pts3c(1,:)));
t = find(pts3c(2,r) == min(pts3c(2,r)));
y = find(pts3c(2,r) == max(pts3c(2,r)));

pts3ct = pts3c(:,q(w));
pts3ct(:,end+1) = pts3c(:,q(e));
pts3ct(:,end+1) = pts3c(:,r(t));
pts3ct(:,end+1) = pts3c(:,r(y));
pts3c = pts3ct;


pgonRange12 = polyshape([pts1c(1,1),pts1c(1,2),pts1c(1,3),pts1c(1,4)],[pts1c(2,1),pts1c(2,2),pts1c(2,3),pts1c(2,4)]);
pgonRange23 = polyshape([pts2c(1,1),pts2c(1,2),pts2c(1,3),pts2c(1,4)],[pts2c(2,1),pts2c(2,2),pts2c(2,3),pts2c(2,4)]);
pgonRange31 = polyshape([pts3c(1,1),pts3c(1,2),pts3c(1,3),pts3c(1,4)],[pts3c(2,1),pts3c(2,2),pts3c(2,3),pts3c(2,4)]);
pgonRange = intersect(pgonRange12,pgonRange23);
pgonRange = intersect(pgonRange,pgonRange31);


%Fulfilling these additional inequalities will give trajectories that allow
%for addition of 1 to 3 more agents (4 to 6 total agents)

%for phi = +pi/2 or -pi/2
%|((k+1)(ri^2 + rj^2)^0.5 - (di^2 + dj^2)^0.5)| > dCT (0.5)

syms dc Rc
% eq = abs(alpha1r)*Rc - abs(alpha1d)*dc == Dmin;
% eqneg = abs(alpha1r)*Rc - abs(alpha1d)*dc == -Dmin;
pts1 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = (k+1)*((alpha1r*Rc/(k+1))^2 + (alpha2r*Rc/(k+1))^2)^0.5 - ((alpha1d*dc)^2 + (alpha2d*dc)^2)^0.5 == CollThresh; 
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eq = (k+1)*((alpha1r*Rc/(k+1))^2 + (alpha2r*Rc/(k+1))^2)^0.5 - ((alpha1d*dc)^2 + (alpha2d*dc)^2)^0.5 == CollThresh;
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV(1);  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = (k+1)*((alpha1r*Rc/(k+1))^2 + (alpha2r*Rc/(k+1))^2)^0.5 - ((alpha1d*dc)^2 + (alpha2d*dc)^2)^0.5 == -CollThresh;
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eqneg = (k+1)*((alpha1r*Rc/(k+1))^2 + (alpha2r*Rc/(k+1))^2)^0.5 - ((alpha1d*dc)^2 + (alpha2d*dc)^2)^0.5 == -CollThresh;
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
            pts1(2,i) = SLV(1);  
        end

    end
end

% eq = abs(alpha2r)*Rc - abs(alpha2d)*dc == Dmin;
% eqneg = abs(alpha2r)*Rc - abs(alpha2d)*dc == -Dmin;
pts2 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = (k+1)*((alpha2r*Rc/(k+1))^2 + (alpha3r*Rc/(k+1))^2)^0.5 - ((alpha2d*dc)^2 + (alpha3d*dc)^2)^0.5 == CollThresh; 
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eq = (k+1)*((alpha2r*Rc/(k+1))^2 + (alpha3r*Rc/(k+1))^2)^0.5 - ((alpha2d*dc)^2 + (alpha3d*dc)^2)^0.5 == CollThresh; 
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV(1);  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = (k+1)*((alpha2r*Rc/(k+1))^2 + (alpha3r*Rc/(k+1))^2)^0.5 - ((alpha2d*dc)^2 + (alpha3d*dc)^2)^0.5 == -CollThresh; 
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eqneg = (k+1)*((alpha2r*Rc/(k+1))^2 + (alpha3r*Rc/(k+1))^2)^0.5 - ((alpha2d*dc)^2 + (alpha3d*dc)^2)^0.5 == -CollThresh; 
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
            pts2(2,i) = SLV(1);  
        end

    end
end

% eq = abs(alpha3r)*Rc - abs(alpha3d)*dc == Dmin;
% eqneg = abs(alpha3r)*Rc - abs(alpha3d)*dc == -Dmin;
pts3 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            eq = (k+1)*((alpha3r*Rc/(k+1))^2 + (alpha1r*Rc/(k+1))^2)^0.5 - ((alpha3d*dc)^2 + (alpha1d*dc)^2)^0.5 == CollThresh; 
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eq = (k+1)*((alpha3r*Rc/(k+1))^2 + (alpha1r*Rc/(k+1))^2)^0.5 - ((alpha3d*dc)^2 + (alpha1d*dc)^2)^0.5 == CollThresh; 
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV(1);  
        end
    else
        if(i==3)
            Rc = CValLim;
            eqneg = (k+1)*((alpha3r*Rc/(k+1))^2 + (alpha1r*Rc/(k+1))^2)^0.5 - ((alpha3d*dc)^2 + (alpha1d*dc)^2)^0.5 == -CollThresh; 
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV(2);
        else
            Rc = -CValLim;
            eqneg = (k+1)*((alpha3r*Rc/(k+1))^2 + (alpha1r*Rc/(k+1))^2)^0.5 - ((alpha3d*dc)^2 + (alpha1d*dc)^2)^0.5 == -CollThresh; 
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
            pts3(2,i) = SLV(1);  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1(1,:) == min(pts1(1,:)));
w = find(pts1(2,q) == max(pts1(2,q)));
e = find(pts1(2,q) == min(pts1(2,q)));
r = find(pts1(1,:) == max(pts1(1,:)));
t = find(pts1(2,r) == min(pts1(2,r)));
y = find(pts1(2,r) == max(pts1(2,r)));

pts1t = pts1(:,q(w));
pts1t(:,end+1) = pts1(:,q(e));
pts1t(:,end+1) = pts1(:,r(t));
pts1t(:,end+1) = pts1(:,r(y));
pts1 = pts1t;

q = find(pts2(1,:) == min(pts2(1,:)));
w = find(pts2(2,q) == max(pts2(2,q)));
e = find(pts2(2,q) == min(pts2(2,q)));
r = find(pts2(1,:) == max(pts2(1,:)));
t = find(pts2(2,r) == min(pts2(2,r)));
y = find(pts2(2,r) == max(pts2(2,r)));

pts2t = pts2(:,q(w));
pts2t(:,end+1) = pts2(:,q(e));
pts2t(:,end+1) = pts2(:,r(t));
pts2t(:,end+1) = pts2(:,r(y));
pts2 = pts2t;

q = find(pts3(1,:) == min(pts3(1,:)));
w = find(pts3(2,q) == max(pts3(2,q)));
e = find(pts3(2,q) == min(pts3(2,q)));
r = find(pts3(1,:) == max(pts3(1,:)));
t = find(pts3(2,r) == min(pts3(2,r)));
y = find(pts3(2,r) == max(pts3(2,r)));

pts3t = pts3(:,q(w));
pts3t(:,end+1) = pts3(:,q(e));
pts3t(:,end+1) = pts3(:,r(t));
pts3t(:,end+1) = pts3(:,r(y));
pts3 = pts3t;


pgonNopi2_1 = polyshape([pts1(1,1),pts1(1,2),pts1(1,3),pts1(1,4)],[pts1(2,1),pts1(2,2),pts1(2,3),pts1(2,4)]);
pgonNopi2_2 = polyshape([pts2(1,1),pts2(1,2),pts2(1,3),pts2(1,4)],[pts2(2,1),pts2(2,2),pts2(2,3),pts2(2,4)]);
pgonNopi2_3 = polyshape([pts3(1,1),pts3(1,2),pts3(1,3),pts3(1,4)],[pts3(2,1),pts3(2,2),pts3(2,3),pts3(2,4)]);
pgonNopi2 = union(pgonNopi2_1,pgonNopi2_2);
pgonNopi2 = union(pgonNopi2,pgonNopi2_3);
pgonpi2 = subtract(pgon,pgonNopi2);


% %For Phi = pi
% %Fulfilling these conditions will allow for the addition of 3 more agents
% %fulfilling phi = pi AND phi = pi/2 will allow for the addition of 9 agents

syms dc Rc

pts1 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            if (sign((alpha1r+alpha2r)*(alpha1d+alpha2d)) == 1)
                eq = (Rc)*(alpha1r+alpha2r) - (alpha1d+alpha2d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha1r+alpha2r) + (alpha1d+alpha2d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
%             pts1(2,i) = S(2);
            pts1(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha1r+alpha2r)*(alpha1d+alpha2d)) == 1)
                eq = (Rc)*(alpha1r+alpha2r) - (alpha1d+alpha2d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha1r+alpha2r) + (alpha1d+alpha2d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts1(1,i) = Rc;
%             pts1(2,i) = S(1);
            pts1(2,i) = SLV; 
        end
    else
        if(i==3)
            Rc = CValLim;
            if (sign((alpha1r+alpha2r)*(alpha1d+alpha2d)) == 1)
                eqneg = (Rc)*(alpha1r+alpha2r) - (alpha1d+alpha2d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha1r+alpha2r) + (alpha1d+alpha2d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
%             pts1(2,i) = S(2);
            pts1(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha1r+alpha2r)*(alpha1d+alpha2d)) == 1)
                eqneg = (Rc)*(alpha1r+alpha2r) - (alpha1d+alpha2d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha1r+alpha2r) + (alpha1d+alpha2d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts1(1,i) = Rc;
%             pts1(2,i) = S(1);
            pts1(2,i) = SLV;
        end

    end
end

% eq = abs(alpha2r)*Rc - abs(alpha2d)*dc == Dmin;
% eqneg = abs(alpha2r)*Rc - abs(alpha2d)*dc == -Dmin;
pts2 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            if (sign((alpha2r+alpha3r)*(alpha2d+alpha3d)) == 1)
                eq = (Rc)*(alpha2r+alpha3r) - (alpha2d+alpha3d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha2r+alpha3r) + (alpha2d+alpha3d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
%             pts2(2,i) = S(2);
            pts2(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha2r+alpha3r)*(alpha2d+alpha3d)) == 1)
                eq = (Rc)*(alpha2r+alpha3r) - (alpha2d+alpha3d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha2r+alpha3r) + (alpha2d+alpha3d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts2(1,i) = Rc;
%             pts2(2,i) = S(1);  
            pts2(2,i) = SLV; 
        end
    else
        if(i==3)
            Rc = CValLim;
            if (sign((alpha2r+alpha3r)*(alpha2d+alpha3d)) == 1)
                eqneg = (Rc)*(alpha2r+alpha3r) - (alpha2d+alpha3d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha2r+alpha3r) + (alpha2d+alpha3d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
%             pts2(2,i) = S(2);
            pts2(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha2r+alpha3r)*(alpha2d+alpha3d)) == 1)
                eqneg = (Rc)*(alpha2r+alpha3r) - (alpha2d+alpha3d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha2r+alpha3r) + (alpha2d+alpha3d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts2(1,i) = Rc;
%             pts2(2,i) = S(1); 
            pts2(2,i) = SLV; 
        end

    end
end

% eq = abs(alpha3r)*Rc - abs(alpha3d)*dc == Dmin;
% eqneg = abs(alpha3r)*Rc - abs(alpha3d)*dc == -Dmin;
pts3 = zeros(2,4);
for i = 1:4
    if(i == 1 || i == 2)
        if(i==1)
            Rc = CValLim;
            if (sign((alpha1r+alpha3r)*(alpha1d+alpha3d)) == 1)
                eq = (Rc)*(alpha1r+alpha3r) - (alpha1d+alpha3d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha1r+alpha3r) + (alpha1d+alpha3d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
%             pts3(2,i) = S(2);
            pts3(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha1r+alpha3r)*(alpha1d+alpha3d)) == 1)
                eq = (Rc)*(alpha1r+alpha3r) - (alpha1d+alpha3d)*dc == CollThresh;
            else
                eq = (Rc)*(alpha1r+alpha3r) + (alpha1d+alpha3d)*dc == CollThresh;
            end
            SLV = solve(eq,dc);
            pts3(1,i) = Rc;
%             pts3(2,i) = S(1); 
            pts3(2,i) = SLV; 
        end
    else
        if(i==3)
            Rc = CValLim;
            if (sign((alpha1r+alpha3r)*(alpha1d+alpha3d)) == 1)
                eqneg = (Rc)*(alpha1r+alpha3r) - (alpha1d+alpha3d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha1r+alpha3r) + (alpha1d+alpha3d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
%             pts3(2,i) = S(2);
            pts3(2,i) = SLV;
        else
            Rc = -CValLim;
            if (sign((alpha1r+alpha3r)*(alpha1d+alpha3d)) == 1)
                eqneg = (Rc)*(alpha1r+alpha3r) - (alpha1d+alpha3d)*dc == -CollThresh;
            else
                eqneg = (Rc)*(alpha1r+alpha3r) + (alpha1d+alpha3d)*dc == -CollThresh;
            end
            SLV = solve(eqneg,dc);
            pts3(1,i) = Rc;
%             pts3(2,i) = S(1);  
            pts3(2,i) = SLV;  
        end

    end
end

pgon = polyshape([-CValLim CValLim CValLim -CValLim],[CValLim,CValLim,-CValLim,-CValLim]);

%Ordering polyshape points
q = find(pts1(1,:) == min(pts1(1,:)));
w = find(pts1(2,q) == max(pts1(2,q)));
e = find(pts1(2,q) == min(pts1(2,q)));
r = find(pts1(1,:) == max(pts1(1,:)));
t = find(pts1(2,r) == min(pts1(2,r)));
y = find(pts1(2,r) == max(pts1(2,r)));

pts1t = pts1(:,q(w));
pts1t(:,end+1) = pts1(:,q(e));
pts1t(:,end+1) = pts1(:,r(t));
pts1t(:,end+1) = pts1(:,r(y));
pts1 = pts1t;

q = find(pts2(1,:) == min(pts2(1,:)));
w = find(pts2(2,q) == max(pts2(2,q)));
e = find(pts2(2,q) == min(pts2(2,q)));
r = find(pts2(1,:) == max(pts2(1,:)));
t = find(pts2(2,r) == min(pts2(2,r)));
y = find(pts2(2,r) == max(pts2(2,r)));

pts2t = pts2(:,q(w));
pts2t(:,end+1) = pts2(:,q(e));
pts2t(:,end+1) = pts2(:,r(t));
pts2t(:,end+1) = pts2(:,r(y));
pts2 = pts2t;

q = find(pts3(1,:) == min(pts3(1,:)));
w = find(pts3(2,q) == max(pts3(2,q)));
e = find(pts3(2,q) == min(pts3(2,q)));
r = find(pts3(1,:) == max(pts3(1,:)));
t = find(pts3(2,r) == min(pts3(2,r)));
y = find(pts3(2,r) == max(pts3(2,r)));

pts3t = pts3(:,q(w));
pts3t(:,end+1) = pts3(:,q(e));
pts3t(:,end+1) = pts3(:,r(t));
pts3t(:,end+1) = pts3(:,r(y));
pts3 = pts3t;


pgonNopi_1 = polyshape([pts1(1,1),pts1(1,2),pts1(1,3),pts1(1,4)],[pts1(2,1),pts1(2,2),pts1(2,3),pts1(2,4)]);
pgonNopi_2 = polyshape([pts2(1,1),pts2(1,2),pts2(1,3),pts2(1,4)],[pts2(2,1),pts2(2,2),pts2(2,3),pts2(2,4)]);
pgonNopi_3 = polyshape([pts3(1,1),pts3(1,2),pts3(1,3),pts3(1,4)],[pts3(2,1),pts3(2,2),pts3(2,3),pts3(2,4)]);
pgonNopi = union(pgonNopi_1,pgonNopi_2);
pgonNopi = union(pgonNopi,pgonNopi_3);
pgonpi = subtract(pgon,pgonNopi);

%Line of cycloids (where r = d)
% eq = (abs(alpha1r)/(k+1))*Rc == abs(alpha1d)*dc;
RcCyclo = 0:CValLim/10:CValLim;
dcCyclo1 = zeros(1,length(RcCyclo));
dcCyclo2 = zeros(1,length(RcCyclo));
dcCyclo3 = zeros(1,length(RcCyclo));

for i = 1:length(RcCyclo)
    dcCyclo1(i) = ((abs(alpha1r)/(k+1))*RcCyclo(i))/abs(alpha1d);
    dcCyclo2(i) = ((abs(alpha2r)/(k+1))*RcCyclo(i))/abs(alpha2d);
    dcCyclo3(i) = ((abs(alpha3r)/(k+1))*RcCyclo(i))/abs(alpha3d);
    
end

%Final area of solutions
FinalPoly = pgonGuard;
FinalPoly = intersect(FinalPoly,pgonMaxDist); %For Max R+r-d
FinalPoly = intersect(FinalPoly,pgonNoColl); %For ensuring no collisions
FinalPoly = intersect(FinalPoly,pgonRange); %For ensuring any 2 agents stay within RangeThresh

FirstQuad = polyshape([0 CValLim CValLim 0],[CValLim,CValLim,0,0]);
FQFinalPoly = intersect(FinalPoly,FirstQuad);
FQFinalPolypi2 = intersect(FQFinalPoly,pgonpi2);
FQFinalPolypi = intersect(FQFinalPoly,pgonpi);
FQFinalPolypipi2 = intersect(FQFinalPolypi,FQFinalPolypi2);
% disp(area(FQFinalPoly)/(CValLim/10))


% %Cleaning out circle cases
% numshapes = 1;
% for i = 1:length(FQFinalPoly.Vertices(:,1))
%     if(isnan(FQFinalPoly.Vertices(i,1)))
%         numshapes = numshapes + 1;
%     end
% end

FQPolyXMax = max(FQFinalPoly.Vertices(:,1));
FQPolyYMax = max(FQFinalPoly.Vertices(:,2));

figure(1); clf(1);
title('Visualization of Solution Space - Geometric Constraints')
plot(FQFinalPoly);
axis equal

figure(2); clf(2);
title('Visualization of Solution Space - Geometric Constraints and Additional Agents')
hold on
plot(FQFinalPoly);
% plot(subtract(FQFinalPolypi2,FQFinalPolypipi2),'FaceColor','g')
plot(subtract(FQFinalPolypi,FQFinalPolypipi2),'FaceColor','r')
plot(FQFinalPolypipi2,'FaceColor','k')
FQPolyMAX = max(FQPolyXMax,FQPolyYMax);
xlim([0 FQPolyMAX])
ylim([0 FQPolyMAX])
legend('Geometric Only', 'Geometric and 3 more agents', 'Geometric and 6 more agents')


%figure(11*ttestno); clf(11*ttestno);
figure(3); clf(3);
title('Visualization of Each Constraint in Solution Space')
plot(pgonRange)
hold on
plot(pgonMaxDist)
plot(pgonGuard)
plot(pgonNoColl)
legend('Range Constraint', 'Maximum Distance Constraint', 'Minimum Distance Constraint', 'Collision Constraint')
% FQPolyXMax = max(FQFinalPoly.Vertices(:,1));
% FQPolyYMax = max(FQFinalPoly.Vertices(:,2));
xlim([0 max(pgonMaxDist.Vertices(:,1))])
ylim([0 max(pgonMaxDist.Vertices(:,1))])
axis equal

% figure(12); clf(12);
% plot(FQFinalPoly)
% hold on
% % plot(RcCyclo,dcCyclo1,'b')
% % plot(RcCyclo,dcCyclo2,'r')
% % plot(RcCyclo,dcCyclo3,'g')
% FQPolyXMax = max(FQFinalPoly.Vertices(:,1));
% FQPolyYMax = max(FQFinalPoly.Vertices(:,2));
% xlim([0 FQPolyXMax])
% ylim([0 FQPolyXMax])

%figure(4); clf(4);
%plot(FQFinalPolypipi2,'FaceColor','k')
%axis equal
%xlim([0 FQPolyMAX])
%ylim([0 FQPolyMAX])


nanc = 0;
for i = 1:length(FQFinalPolypipi2.Vertices)
    if isnan(FQFinalPolypipi2.Vertices(i,1))
        nanc = nanc + 1;
    end
end

numshape = nanc + 1;

scale3AFact = (max(max(FQFinalPoly.Vertices)));
scaled3APoly = polyshape(FQFinalPoly.Vertices/(scale3AFact));

scaleFact = (max(max(FQFinalPolypipi2.Vertices)));
scaledPoly = polyshape(FQFinalPolypipi2.Vertices/(scaleFact));
% scaledPoly = polyshape(FQFinalPolypipi2.Vertices/100);

figure(4); clf(4);
title('Scaled Down Solution Space - Geometric Constraints Only')
hold on
plot(scaled3APoly)
axis equal

figure(5); clf(5);
plot(scaledPoly)
hold on
title('Scaled Down Solution Space - Geometric Constraints and 6 Extra Agents')

axis equal
xlim([0 max(max(scaledPoly.Vertices))])
ylim([0 max(max(scaledPoly.Vertices))])

cntrdAng = zeros(4,numshape);
newshapepos = [1; find(isnan(scaledPoly.Vertices(:,1)))+1];

% currshape = 1;
for currshape = 1:numshape
    
    sumx = 0;
    sumy = 0;

    st = newshapepos(currshape);
    if currshape == numshape
        en = size(scaledPoly.Vertices,1);
    else
        en = newshapepos(currshape+1) - 2;
    end

    for j = st:en
        sumx = sumx + scaledPoly.Vertices(j,1);
        sumy = sumy + scaledPoly.Vertices(j,2);
    end
    cntrdAng(1,currshape) = sumx/(en-st+1);
    cntrdAng(2,currshape) = sumy/(en-st+1);
    cntrdAng(3,currshape) = atan2(cntrdAng(2,currshape),cntrdAng(1,currshape));
    cntrdAng(4,currshape) = abs(cntrdAng(3,currshape)-pi/4);

end

LCRes = find(cntrdAng(4,:)==min(cntrdAng(4,:)));
RelevST = newshapepos(LCRes);
if (LCRes == numshape)
    RelevEN = size(scaledPoly.Vertices,1);
else
    RelevEN = newshapepos(LCRes+1) - 2;
end


relevx = scaledPoly.Vertices(RelevST:RelevEN,1);
relevy = scaledPoly.Vertices(RelevST:RelevEN,2);
scaledRPoly = polyshape(relevx,relevy);

RArea = area(scaledRPoly);

polytest = scaledRPoly;
% polytest = polyshape([1 2 1.5 -1],[0 0.5 1.75 1]);


dimno = 2;
c_mat = zeros(length(polytest.Vertices(:,1)),1);
b_mat = zeros(length(polytest.Vertices(:,1)),1);
A_mat = zeros(length(polytest.Vertices(:,1)),dimno);
pts = polytest.Vertices;
m_mat = zeros(length(polytest.Vertices(:,1)),1);
xm = c_mat;

pt1 = polytest.Vertices(1,:);
pt2 = polytest.Vertices(2,:);
pt3 = polytest.Vertices(3,:);

polycentroid = [sum(polytest.Vertices(:,1))/length(polytest.Vertices(:,1)) sum(polytest.Vertices(:,2))/length(polytest.Vertices(:,2))];

for i = 1:size(pts,1)
    if(i == size(pts,1))
        m_mat(i) = (pts(1,2)-pts(i,2))/(pts(1,1)-pts(i,1));
    else
        m_mat(i) = (pts(i+1,2)-pts(i,2))/(pts(i+1,1)-pts(i,1));
    end

    if(abs(m_mat(i))==Inf)
        xm(i) = 1;
    else
        xm(i) = 0;
    end
    
    if(xm(i)==0)
        c_mat(i) = -m_mat(i)*pts(i,1) + pts(i,2);

        if(-m_mat(i)*polycentroid(1) + polycentroid(2) < c_mat(i))
            b_mat(i) = c_mat(i);
            A_mat(i,:) = [-m_mat(i) 1];
        else
            b_mat(i) = -c_mat(i);
            A_mat(i,:) = -[-m_mat(i) 1];
        end
    else
        c_mat(i) = -pts(i,1) + (1/m_mat(i))*pts(i,2);

        if(-polycentroid(1) + (1/m_mat(i))*polycentroid(2) < c_mat(i))
            b_mat(i) = c_mat(i);
            A_mat(i,:) = [-1 (1/m_mat(i))];
        else
            b_mat(i) = -c_mat(i);
            A_mat(i,:) = -[-1 (1/m_mat(i))];
        end
    end
    
    
end

A = double(A_mat);
b = double(b_mat);

m = size(pts,1);
n = dimno;
cvx_begin
    variable Bel(n,n) symmetric
    variable del(n)
    maximize( det_rootn(Bel))
    subject to
        for i=1:m
            norm(Bel*A(i,:)',2) + A(i,:)*del <= b(i);
        end
cvx_end

ielli = Bel*[cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))] + del*ones(1,100);
ielli2 = Bel/2*[cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))] + del*ones(1,100);
ielli3 = Bel/3*[cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))] + del*ones(1,100);

ielli2scaled = ielli2*scaleFact;
ielli3scaled = ielli3*scaleFact;
elli_dist = zeros(1,length(ielli));

% minaxcirc = [min(eig(Bel))*cos(linspace(0,2*pi,100)); min(eig(Bel))*sin(linspace(0,2*pi,100))] + del*ones(1,100);
% minaxcircscaled = minaxcirc;

for i = 1:length(ielli)
    elli_dist(i) = norm([del(1) del(2)]-[ielli(1,i) ielli(2,i)]);
end

ellimax = find(elli_dist==max(elli_dist));
ellimin = find(elli_dist==min(elli_dist));

minax = norm([del(1) del(2)]-[ielli(1,ellimin) ielli(2,ellimin)]);
majax = norm([del(1) del(2)]-[ielli(1,ellimax) ielli(2,ellimax)]);

%figure(888); clf(888);
%plot([1:1:length(elli_dist)],elli_dist);

figure(6); clf(6);
plot(polytest)
axis equal
hold on
plot(polycentroid(1),polycentroid(2),'ok')
plot(ielli(1,:),ielli(2,:),'r')
plot(del(1),del(2),'or')
plot([del(1) ielli(1,ellimax)],[del(2) ielli(2,ellimax)],'r')
plot([del(1) ielli(1,ellimin)],[del(2) ielli(2,ellimin)],'r')

plot(ielli2(1,:),ielli2(2,:),'b')
plot(ielli3(1,:),ielli3(2,:),'g')

dscaled = del*scaleFact;

C01 = dscaled';
syms x10 x20 y10 y20
PhiA = iPhiA; PhiB = iPhiB;
eqx1 = C01(1) == sqrt((G1*x10 + G2*x20)^2 + (G1*y10 + G2*y20)^2);
eqx2 = C01(2) == sqrt((G1p*x10 + G2p*x20)^2 + (G1p*y10 + G2p*y20)^2);
eqx3 = atan2((G1Phi*y10 + G2Phi*y20),(G1Phi*x10 + G2Phi*x20)) == PhiA;
eqx4 = atan2((G1Phip*y10 + G2Phip*y20),(G1Phip*x10 + G2Phip*x20)) == PhiB;

HardSolve = solve([eqx1,eqx2,eqx3,eqx4],[x10, x20, y10, y20]);
x0 = [round(HardSolve.x10,3);round(HardSolve.y10,3);round(HardSolve.x20,3);round(HardSolve.y20,3);0;0];
x0_Orig = x0;

M = -(kron((G*L),eye(2,2))-kron((B*L),S));
V=null(M');
M1 = [V(:,1)'*(kron(ones(3,1),eye(2,2))); V(:,2)'*(kron(ones(3,1),eye(2,2)))];
M2 = [V(:,1)'*x0; V(:,2)'*x0];
xInf = M1\M2;

x0 = x0 - repmat(xInf,3,1);
x0 = round(x0,3);

M = -(kron((G*L),eye(2,2))-kron((B*L),S));
V=null(M');
M1 = [V(:,1)'*(kron(ones(3,1),eye(2,2))); V(:,2)'*(kron(ones(3,1),eye(2,2)))];
M2Shift = [V(:,1)'*x0; V(:,2)'*x0];
xInfShift = M1\M2Shift;

x0xarr = [x0(1); x0(3); x0(5);];
x0yarr = [x0(2); x0(4); x0(6);];


PolyDiscr = 100;
PolyPrimeX = linspace(min(scaledRPoly.Vertices(:,1)),max(scaledRPoly.Vertices(:,1)),PolyDiscr);
PolyPrimeY = linspace(min(scaledRPoly.Vertices(:,2)),max(scaledRPoly.Vertices(:,2)),PolyDiscr);
PolyPrimeArr = zeros(2,PolyDiscr^2);

for m = 1:length(PolyPrimeX)
    for n = 1:length(PolyPrimeY)
        PolyPrimeArr(:,(m-1)*PolyDiscr + n) = [PolyPrimeX(m); PolyPrimeY(n)];
    end
end

PolyPrimeFArr = PolyPrimeArr;
FArrCount = 1;

for l = 1:length(PolyPrimeArr)
    if(inpolygon(PolyPrimeArr(1,l),PolyPrimeArr(2,l),scaledRPoly.Vertices(:,1),scaledRPoly.Vertices(:,2))==1)
        PolyPrimeFArr(:,FArrCount) = PolyPrimeArr(:,l);
        FArrCount = FArrCount + 1;
    end
end

PolyPrimeFFArr = PolyPrimeFArr(:,1:FArrCount-1);

%let PhiA, PhiB vary 
PhiMin = -pi/16;
PhiMax = pi/16;

%Find Cshift
mind1 = abs(abs(alpha1r)*dscaled(1) - abs(alpha1d)*dscaled(2));
mind2 = abs(abs(alpha2r)*dscaled(1) - abs(alpha2d)*dscaled(2));
mind3 = abs(abs(alpha3r)*dscaled(1) - abs(alpha3d)*dscaled(2));

maxd1 = abs(abs(alpha1r)*dscaled(1) + abs(alpha1d)*dscaled(2));
maxd2 = abs(abs(alpha2r)*dscaled(1) + abs(alpha2d)*dscaled(2));
maxd3 = abs(abs(alpha3r)*dscaled(1) + abs(alpha3d)*dscaled(2));

mindd1 = mind1 - Dmin;
mindd2 = mind2 - Dmin;
mindd3 = mind3 - Dmin;
maxdd1 = Dmax - maxd1;
maxdd2 = Dmax - maxd2;
maxdd3 = Dmax - maxd3;

ShiftCompArr = [mindd1 mindd2 mindd3 maxdd1 maxdd2 maxdd3];
Cshift = min(ShiftCompArr);

