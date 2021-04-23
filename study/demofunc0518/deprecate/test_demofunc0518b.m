myclear;
%
funcPrm.numerPrm.c10 = 0.0;
funcPrm.numerPrm.numTerms1 = 2;
funcPrm.numerPrm.c1(1,:) = [ 1.0, 3.0, 0, 0.0 ];
funcPrm.numerPrm.c1(2,:) = [ 2.0, 1.0, 1, 0.0 ];
%
funcPrm.numerPrm.c20 = 0.0;
funcPrm.numerPrm.numTerms2 = 0;
funcPrm.numerPrm.c2 = [ 0.0, 1.0, 1, pi ];
%
funcPrm.denomPrm.c10 = 1.0;
funcPrm.denomPrm.numTerms1 = 0;
funcPrm.denomPrm.c1 = [ 1.0, 2.0, 0, 0.0 ];
%
funcPrm.denomPrm.c20 = 1.0;
funcPrm.denomPrm.numTerms2 = 0;
funcPrm.denomPrm.c2 = [ 0.0, 2.0, 0, 0.0 ];
%
%
funchF = @(x) demofunc0518b_eval(x,funcPrm);
funchF1XY = @(x,y) funchF([x;y])(1,:);
funchF2XY = @(x,y) funchF([x;y])(2,:);
funchFNXY = @(x,y) sqrt(sum(funchF([x;y]).^2,1));
%
funchVizXY = @(x,y) sqrt(funchFNXY(x,y));
%
contourfunch(funchVizXY,2*[-1,1,-1,1]);
axis equal;
cMap = 0.6 + (0.4*jet);
cMap(1,:) *= 0.4;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
colormap(cMap);
grid on;
