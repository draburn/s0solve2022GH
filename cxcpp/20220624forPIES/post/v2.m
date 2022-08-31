numFigs = 0;
jAbsMax = max(max(abs(matJ)));
%
matXIndex = repmat( (1:sizeX), [sizeF,1] );
matFIndex = repmat( (1:sizeF)', [1,sizeX] );
%
if (0)
sz = 30000;
alpha0 = linspace( 0.0, 1.0, sz )';
alphaC = alpha0.*(1.0-alpha0);
alphaC = alphaC.^2;
alphaC /= max(alphaC);
alphaL = (1.0-alphaC).*(alpha0<0.5);
alphaR = (1.0-alphaC).*(alpha0>0.5);
cmapL = zeros(sz,3);
cmapC = 0.5-0.4*jet(sz);
foo = cmapC(:,1);
cmapC(:,1) = cmapC(:,3);
cmapC(:,3) = foo;
cmapC(:,2) *= 0.5;
cmapR = ones(sz,3);
cmap1 = (cmapL.*alphaL) + (cmapC.*alphaC) + (cmapR.*alphaR); % Automatic broadcast.
cmap1 = 1.0 - cmap1;
%cmap1 = cmap1./( 0.6*cmap1(:,1) + cmap1(:,2) + 0.1*cmap1(:,3) );  % Automatic broadcast.
cmap1 = (1.0 + (cmap1-1.0).*(cmap1<1.0)).*(cmap1>0.0);
cmapPos = cmap1;
%
numFigs++;figure(numFigs);
imagesc( matXIndex );
title( "matXIndex" );
colormap(cmapPos);
%
numFigs++;figure(numFigs);
imagesc( matFIndex );
title( "matFIndex" );
colormap(cmapPos);
%
numFigs++;figure(numFigs);
imagesc( asinh(abs(matJ)/jAbsMax) );
title( "asinh(abs(matJ)/jAbsMax)" );
colormap(cmapPos);
%
numFigs++;figure(numFigs);
imagesc( asinh(1.0e4*abs(matJ)/jAbsMax)/1.0e4 );
title( "asinh(1.0e4*abs(matJ)/jAbsMax)/1.0e4" );
colormap(cmapPos);
endif
%
sz = 30000;
alpha0 = linspace( 0.0, 1.0, sz )';
alphaC = alpha0.*(1.0-alpha0);
alphaC = alphaC.^2;
alphaC /= max(alphaC);
alphaL = (1.0-alphaC).*(alpha0<0.5);
alphaR = (1.0-alphaC).*(alpha0>0.5);
cmapL = [ linspace(0.0,0.0,sz)', linspace(0.0,1.0,sz)'.^0.3, linspace(1.0,0.0,sz)' ];
cmapC = ones(sz,3);
cmapR = [ linspace(0.0,1.0,sz)', linspace(1.0,0.0,sz)'.^0.3, linspace(0.0,0.0,sz)' ];
cmap1 = (cmapL.*alphaL) + (cmapC.*alphaC) + (cmapR.*alphaR); % Automatic broadcast.
%cmap1 = 1.0 - cmap1;
%cmap1 = cmap1./( 0.6*cmap1(:,1) + cmap1(:,2) + 0.1*cmap1(:,3) );  % Automatic broadcast.
cmap1 = (1.0 + (cmap1-1.0).*(cmap1<1.0)).*(cmap1>0.0);
cmapCent = cmap1;
%
numFigs++;figure(numFigs);
image( sz*(1+matXIndex)/(1.0+sizeX) );
title( "((matXIndex-1.0)/(sizeX-1.0))-0.5" );
colormap(cmapCent);
%
numFigs++;figure(numFigs);
image( sz*(1+matFIndex)/(1.0+sizeF) );
title( "((matFIndex-1.0)/(sizeF-1.0))-0.5" );
colormap(cmapCent);
%
numFigs++;figure(numFigs);
image( sz*0.5*(1.0+asinh(matJ/jAbsMax)/asinh(1.0)) );
title( "asinh(matJ/jAbsMax)/asinh(1.0)" );
colormap(cmapCent);
%
numFigs++;figure(numFigs);
image( sz*0.5*(1.0+asinh(1.0e6*matJ/jAbsMax)/asinh(1.0e6)) );
title( "asinh(1.0e3*matJ/jAbsMax)/asinh(1.0e3)" );
colormap(cmapCent);
%
numFigs++;figure(numFigs);
image( sz*0.5*(1.0+asinh(1.0e6*matJ/jAbsMax)/asinh(1.0e6)) );
title( "asinh(1.0e6*matJ/jAbsMax)/asinh(1.0e6)" );
colormap(cmapCent);
%
%
condJ = cond(matJ)
vecS = svd(matJ);
%
numFigs++;figure(numFigs);
semilogy( vecS, 'o-', 'markersize', 20, 'linewidth', 2 );
grid on;
title( 'singular values of matJ' );

return
%
matH = matJ'*matJ;
[ matPsi, matLambda ] = eig( matH );
vecLambda = diag(matLambda);
%
numFigs++;figure(numFigs);
semilogy( vecLambda, 'o-', 'markersize', 20, 'linewidth', 2 );
grid on;
title( 'eigenvalues of Hessian' );
