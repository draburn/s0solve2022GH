% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
%
sizeF = 1;
sizeX = 1000;
numElemPerCol = 0;
numAddtlElemPerRow = 5;
c0 = 0.00;
csx = 0.0;
csf = 0.0;
%
matJ0 = zeros(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	%matJ0(m,n) = randn();
	matJ0(m,n) = 1.0;
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	%matJ0(m,n) = randn();
	matJ0(m,n) = 1.0;
	actualElemUsed(t) = n;
endfor
endfor
matJ0 += c0 * randn(sizeF,sizeX);
assert( isrealarray(matJ0,[sizeF,sizeX]) );
actualElemUsed

matSX = diag(exp(csx*randn(sizeX,1)));
matSF = diag(exp(csf*randn(sizeF,1)));
matJ = matSF*matJ0/matSX;


%setprngstates(3568384);
%setprngstates(74504256);
%setprngstates(26048640);
setprngstates(38115184);
sizeU0 = 55;

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(matU0);

sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;

matJL2_1 = (matW*(matV'))/( matV*(matV') + 1.0*sqrt(eps)*diag(diag(matV*(matV'))) );
matJL2_2 = (matW*(matV'))/( matV*(matV') + 2.0*sqrt(eps)*diag(diag(matV*(matV'))) );
matJL2 = 2.0*matJL2_1 - matJL2_2;
%
if (1)
%error( "This works better iteratively; see calcSparseMatrixEstimate.m." ); No that's not quite the same!
matJL2E = zeros(sizeF,sizeX);
for m=1:sizeF
	[ foo, maxyElem ] = sort(-abs(matJL2(m,:)));
	elemSelector = maxyElem(1:sizeV);
	%elemSelector = maxyElem(1:ceil(sizeV/2.0));
	%elemSelector = maxyElem(1:5);
	matVSel = matV(elemSelector,:);
	matJL2E(elemSelector) = matW(m,:)*(matVSel')/( matVSel*(matVSel') );
endfor
endif

matJEst = calcSparseMatrixEstimate( matV, matW );

%
m = 1;
numFigs++; figure(numFigs);
plot( matJ(m,:), 'o-', 'linewidth', 5, matJL2E(m,:), 's-', 'markersize', 10, 'linewidth', 1, matJEst(m,:), 'x-', 'linewidth', 2 );
grid on;
return
numFigs++; figure(numFigs);
plot( matChi(m,:), '^-', 'linewidth', 2 );
grid on;
numFigs++; figure(numFigs);
plot( matRho(m,:), '^-', 'linewidth', 2 );
grid on;




return;
