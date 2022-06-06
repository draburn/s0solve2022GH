% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates();
tic();
%
sizeF = 100;
sizeX = 100;
numElemPerCol = 0;
numAddtlElemPerRow = 5;
c0 = 0.01;
csx = 0.0;
csf = 0.0;
%
matJ0 = zeros(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
	%matJ0(m,n) = 1.0;
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
	%matJ0(m,n) = 1.0;
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
%setprngstates(38115184);
sizeU0 = 50;
%sizeU0 = 80;

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(matU0);

sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;

matJL2 = zeros(sizeF,sizeX);
if (0)
matJL2_1 = (matW*(matV'))/( matV*(matV') + 1.0*sqrt(eps)*diag(diag(matV*(matV'))) );
matJL2_2 = (matW*(matV'))/( matV*(matV') + 2.0*sqrt(eps)*diag(diag(matV*(matV'))) );
matJL2 = 2.0*matJL2_1 - matJL2_2;
endif
%
if (0)
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

toc();
t0=time();
matJEst_hacky = calcSparseMatrixEstimate( matV, matW );
msg( __FILE__, __LINE__, sprintf( "matJEst_hacky() took %0.3gms", 1000.0*(time()-t0)) );
t0=time();
matJEst_basic = calcSparseMatEst_basic( matV, matW );
msg( __FILE__, __LINE__, sprintf( "matJEst_basic() took %0.3gms", 1000.0*(time()-t0)) );

%
m = 1;
numFigs++; figure(numFigs);
plot( ...
  matJ(m,:), 'o-', 'linewidth', 6, 'markersize', 5, ...
  matJL2(m,:), 's-', 'linewidth', 1, 'markersize', 20, ...
  matJEst_hacky(m,:), 'x-', 'linewidth', 4, 'markersize', 10, ...
  matJEst_basic(m,:), '^-', 'linewidth', 2, 'markersize', 15 );
grid on;

numFigs++; figure(numFigs);
imagesc(abs(matJ));
grid on;
axis equal;
numFigs++; figure(numFigs);
imagesc(abs(matJEst_basic));
grid on;
axis equal;
numFigs++; figure(numFigs);
imagesc(matJEst_basic-matJ);
grid on;
axis equal;
max(max(abs(matJEst_basic-matJ)))

return
