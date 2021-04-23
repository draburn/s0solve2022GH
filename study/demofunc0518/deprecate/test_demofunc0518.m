myclear;
t1 = 0;
t2 = 0;
c1 = [];
c2 = [];
t1++; c1(t1,:) = [ 1.0, 3.0, 0, 0.0 ];
t2++; c2(t2,:) = [ 1.0, 3.0, 0, 0.0 ];
t1++; c1(t1,:) = [ 2.0, 2.0, 3, 0.0 ];
%t2++; c1(t2,:) = [ 1.0, 1.0, 3, 0.0 ];
%
funcPrm.numTerms1 = t1;
funcPrm.numTerms2 = t2;
funcPrm.c10 = 0.0;
funcPrm.c20 = 0.0;
funcPrm.c1 = c1;
funcPrm.c2 = c2;
%
funchF = @(x) demofunc0518_eval(x,funcPrm);
funchF1XY = @(x,y) funchF([x;y])(1);
funchF2XY = @(x,y) funchF([x;y])(2);
funchFNXY = @(x,y) sqrt(sum(funchF([x;y]).^2,1));
%
funchVizXY = @(x,y) sqrt(funchFNXY(x,y));
%
contourfunch(funchVizXY,[-1,1,-1,1])
axis equal;
cMap = 0.6 + (0.4*jet);
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
colormap(cMap);
grid on;
