% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(71021024);
%
sizeF = 1;
sizeX = 100;
numElemPerCol = 0;
numAddtlElemPerRow = 6;
c0 = 0.0;
csx = 0.0;
csf = 0.0;
%
matJ0 = zeros(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
matJ0 += c0 * randn(sizeF,sizeX);
assert( isrealarray(matJ0,[sizeF,sizeX]) );

matSX = diag(exp(csx*randn(sizeX,1)));
matSF = diag(exp(csf*randn(sizeF,1)));
matJ = matSF*matJ0/matSX;



%setprngstates(96551984); % Fail
setprngstates(75868048); % Succ
sizeU0 = 90;

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(randn(sizeX,sizeU0));
if (0)
	matU3 = zeros(sizeX,3);
	matU3(1:3:end,1) = 1.0;
	matU3(2:3:end,2) = 1.0;
	matU3(3:3:end,3) = 1.0;
	matU5 = zeros(sizeX,5);
	matU5(1:5:end,1) = 1.0;
	matU5(2:5:end,2) = 1.0;
	matU5(3:5:end,3) = 1.0;
	matU5(4:5:end,4) = 1.0;
	matU5(5:5:end,5) = 1.0;
	matV = [ matV, matU3, matU5 ];
endif

%matV = matV(:,1:10);
%matV = matV(:,11:20);

sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;
%
%
%
matJCollectiveEst = zeros(sizeF,sizeX);
matResCollective = zeros(sizeF,sizeV);
for m=1:sizeF
	funchOmega = @(x)( 1E16*norm(matV'*x-matW(m,:)') + norm(x,1) );
	vecX0 = zeros(sizeX,1);
	fminunc_options = optimset( "TolFun", 1e-20 );
	[ vecXF, fminunc_fval, fminunc_info, fminunc_output, fminunc_grad, fminunc_hess ] = fminunc(funchOmega,vecX0,fminunc_options);
	fminunc_fval
	fminunc_output
	fminunc_info
	funchOmega(vecXF)
	matJEst(m,:) = vecXF;
endfor
matResCollective = matW - matJEst*matV;
%matJCollectiveEst
%matJ

if (1)
	numFigs++; figure(numFigs);
	plot( matJ, 'o-', 'linewidth', 3, matJEst, 'x-' );
	grid on;
	numFigs++; figure(numFigs);
	plot( matResCollective, 'o-' );
	grid on;
endif
