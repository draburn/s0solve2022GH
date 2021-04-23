myclear;
setprngstates(0);
%
sizeN = 2;
%
if (0)
	matA = eye(sizeN,sizeN);
	matB = eye(sizeN,sizeN);
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
elseif (1)
	%matA = [ 1.0, 0.01; 1.0, 0.5 ];
	%matB = [ 1.0, 0.0; -0.1, 1.0 ];
	matA = [ 1.0, 0.01; 1.0, 0.1 ];
	matB = [ 1.0, 0.0; -0.1, 1.0 ];
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
elseif (1)
	matA = eye(sizeN,sizeN) + 0.5*randn(sizeN,sizeN);
	matB = eye(sizeN,sizeN) + 0.5*randn(sizeN,sizeN);
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
else
	matA = randn(sizeN,sizeN);
	matB = randn(sizeN,sizeN);
	vecX0 = randn(sizeN,1);
	indexNNonLin = ceil((10*eps) + ((sizeN-(20*eps))*rand) );
end
%
c1 = 1.0;
cx = 0.001;
c0 = (c1^2)*(1.0+cx)/4.0;
funchF = @(x)(blm0429_func( x, matA, matB, vecX0, indexNNonLin, c1, c0 ));
