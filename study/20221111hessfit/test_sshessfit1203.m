clear;
mydefs;
setprngstates(0);
%
sizeX = 10
gNoiseLevel = 0.0
fNoiseLevel = 0.0
numPts = 5
%
f0 = 0.0;
vecX0 = randn(sizeX,1);
sizeL = ceil(sizeX*sqrt(sizeX));
matA = sparse(diag(randn(sizeX,1))) + sparse( ...
  1 + floor( (sizeX-eps) * rand(sizeL,1) ), ...
  1 + floor( (sizeX-eps) * rand(sizeL,1) ), ...
  randn(sizeL,1), ...
  sizeX, ...
  sizeX );
sizeTemp = 1;
matTemp = randn(sizeX,sizeTemp);
matA /= sqrt( sizeX * sum(sumsq(matA'*matTemp,1) ./ sumsq(matTemp,1),2) / sizeTemp );
%
funchD = @(x)( x - vecX0 );
funchHDOfD = @(d)( matA*((matA')*d) );
funchFSmoothOfD = @(d)( f0 + 0.5*sum( d .* funchHDOfD(d), 1 ) );
funchGSmoothOfD = @(d)( funchHDOfD(d) );
funchFSmooth = @(x)( funchFSmoothOfD(funchD(x)) );
funchGSmooth = @(x)( funchGSmoothOfD(funchD(x)) );
funchF = @(x)( funchFSmooth(x) + fNoiseLevel*randn([1,size(x,2)]) );
funchG = @(x)( funchGSmooth(x) + gNoiseLevel*randn(size(x)) );
%
matX = randn(sizeX,numPts);
rvecF = funchF(matX);
matG = funchG(matX);
%
switch (10)
case 0
	vecA = matX(:,1);
	matV = eye(sizeX);
	fSS = funchF(vecA);
	vecGSS = funchG(vecA);
	msg( __FILE__, __LINE__, "Calculating full Hessian..." ); tic();
	matHSS = matA * (matA');
	msgnnl( __FILE__, __LINE__, "Finished calculating full hessian. " ); toc();
case 9
	prm = [];
	msg( __FILE__, __LINE__, "Calling sshessfit1202()..." ); tic();
	[ vecA, matV, fSS, vecGSS, matHSS ] = sshessfit1202( sizeX, numPts, matX, rvecF, matG, prm );
	msgnnl( __FILE__, __LINE__, "Finished sshessfit1202(). " ); toc();
case 10
	prm = [];
	msg( __FILE__, __LINE__, "Calling sshessfit1203()..." ); tic();
	[ vecA, matV, fSS, vecGSS, matHSS ] = sshessfit1203( sizeX, numPts, matX, rvecF, matG, prm );
	msgnnl( __FILE__, __LINE__, "Finished sshessfit1203(). " ); toc();
case 20
	prm = [];
	msg( __FILE__, __LINE__, "Calling hessfit()..." ); tic();
	[ fSS, vecGSS, matHSS ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm );
	msgnnl( __FILE__, __LINE__, "Finished hessfit(). " ); toc();
	vecA = zeros(sizeX,1);
	matV = eye(sizeX);
case 30
	[ foo, pt0 ] = min(rvecF)
	prm = [];
	msg( __FILE__, __LINE__, "Calling sshessfit1115_orth()..." ); tic();
	[ vecA, matV, fSS, vecGSS, matHSS ] = sshessfit1115_orth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
	msgnnl( __FILE__, __LINE__, "Finished sshessfit1115_orth(). " ); toc();
endswitch
sizeK = size(matV,2);
vecLambda = eig(matHSS);
msg( __FILE__, __LINE__, sprintf( "Eigen value range: %11.3e ~ %11.3e", min(vecLambda), max(vecLambda) ) );
if (min(vecLambda)<=eps)
	matHSS += ( abs(min(vecLambda)) + sqrt(eps)*max(abs(vecLambda)) ) * eye(sizeK);
endif

matRSS = chol( matHSS );
vecDeltaSS = matRSS \ ( matRSS' \ (-vecGSS) );
vecXNext = vecA + (matV * vecDeltaSS);
%
msg( __FILE__, __LINE__, "Results..." );
[ foo, ptAIdeal ] = min(sumsq(matX-vecX0,1));
vecAIdeal = matX(:,ptAIdeal);
msg( __FILE__, __LINE__, sprintf( "  ||x-x_0||: (%10.3e) %10.3e -> %10.3e", norm(vecAIdeal-vecX0), norm(vecA-vecX0), norm(vecXNext-vecX0) ) );
msg( __FILE__, __LINE__, sprintf( "          f: (%10.3e) %10.3e -> %10.3e", funchFSmooth(vecAIdeal), funchFSmooth(vecA), funchFSmooth(vecXNext) ) );
msg( __FILE__, __LINE__, sprintf( "      ||g||: (%10.3e) %10.3e -> %10.3e", norm(funchGSmooth(vecAIdeal)), norm(funchGSmooth(vecA)), norm(funchGSmooth(vecXNext)) ) );
msg( __FILE__, __LINE__, sprintf( "  ||V^T*g||: (%10.3e) %10.3e -> %10.3e", norm(matV'*funchGSmooth(vecAIdeal)), norm(matV'*funchGSmooth(vecA)), norm(matV'*funchGSmooth(vecXNext)) ) );

return;
