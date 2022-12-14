clear;
mydefs;
setprngstates(0);
%
sizeX = 1000
numPts = 21
numSubPtsPerSuperPt = 100
gNoiseLevel = 0.01
fNoiseLevel = 0.01
%
f0 = 0.0;
vecX0 = randn(sizeX,1);
% Or, just use sprandn()?
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
if ( 1 == numSubPtsPerSuperPt )
	matX = randn(sizeX,numPts);
	rvecF = funchF(matX);
	matG = funchG(matX);
else
	avgOverSubPts = @(z)( sum(z,2)/numSubPtsPerSuperPt );
	matX = zeros(sizeX,numPts);
	rvecF = zeros(1,numPts);
	matG = zeros(sizeX,numPts);
	for p = 1 : numPts
		vecXTemp = randn(sizeX,1);
		matXTemp = vecXTemp + (0.1*randn(sizeX,numSubPtsPerSuperPt));
		rvecFTemp = funchF(matXTemp);
		matGTemp = funchG(matXTemp);
		matX(:,p) = avgOverSubPts( matXTemp );
		matG(:,p) = avgOverSubPts( matGTemp );
		%rvecF(p) = sum(rvecFTemp,2)/numSubPtsPerSuperPt;
		rvecF(p) = avgOverSubPts(rvecFTemp) - 0.5*( avgOverSubPts(sum(matXTemp.*matGTemp,1)) - matX(:,p)'*matG(:,p) );
		%[ avgOverSubPts(rvecFTemp), rvecF(p) ]
		clear vecXTemp;
		clear matXTemp;
		clear rvecFTemp;
		clear matGTemp;
	endfor
endif
%
switch (20)
case 0
	vecA = matX(:,1);
	matV = eye(sizeX);
	f = funchFSmooth(vecA);
	vecG = funchGSmooth(vecA);
	msg( __FILE__, __LINE__, "Calculating full Hessian..." ); tic();
	matH = matA * (matA');
	msgnnl( __FILE__, __LINE__, "Finished calculating full hessian. " ); toc();
case 10
	prm = [];
	msg( __FILE__, __LINE__, "Calling hessfit()..." ); tic();
	[ f, vecG, matH ] = hessfit( matX, rvecF, matG, prm );
	msgnnl( __FILE__, __LINE__, "Finished hessfit(). " ); toc();
	vecA = zeros(sizeX,1);
	matV = eye(sizeX);
case 20
	prm = [];
	msg( __FILE__, __LINE__, "Calling sshessfit()..." ); tic();
	[ vecA, matV, f, vecG, matH ] = sshessfit( matX, rvecF, matG, prm );
	msgnnl( __FILE__, __LINE__, "Finished sshessfit(). " ); toc();
endswitch
sizeK = size(matV,2);
vecLambda = eig(matH);
if (min(vecLambda)<=eps)
	matH += ( abs(min(vecLambda)) + sqrt(eps)*max(abs(vecLambda)) ) * eye(sizeK);
endif

matR = chol( matH );
vecY = matR \ ( matR' \ (-vecG) );
vecXNext = vecA + (matV * vecY);
fModelNext = f + vecG'*vecY + 0.5*(vecY'*matH*vecY);
vecGModelNext = matV*(vecG + matH*vecY);
%
msg( __FILE__, __LINE__, "General results..." );
[ foo, ptAIdeal ] = min(sumsq(matX-vecX0,1));
vecAIdeal = matX(:,ptAIdeal);
msg( __FILE__, __LINE__, sprintf( "  sizeK = %d / %d", sizeK, sizeX ) );
msg( __FILE__, __LINE__, sprintf( "  Eigen value range: %11.3e ~ %11.3e", min(vecLambda), max(vecLambda) ) );
msg( __FILE__, __LINE__, "Step..." );
msg( __FILE__, __LINE__, sprintf( "  ||x-x_0||: (%10.3e) %10.3e -> %10.3e (%10.3e)", norm(vecAIdeal-vecX0), norm(vecA-vecX0), norm(vecXNext-vecX0), 0.0 ) );
msg( __FILE__, __LINE__, sprintf( "          f: (%10.3e) %10.3e -> %10.3e (%10.3e)", funchFSmooth(vecAIdeal), funchFSmooth(vecA), funchFSmooth(vecXNext), fModelNext ) );
msg( __FILE__, __LINE__, sprintf( "      ||g||: (%10.3e) %10.3e -> %10.3e (%10.3e)", norm(funchGSmooth(vecAIdeal)), norm(funchGSmooth(vecA)), norm(funchGSmooth(vecXNext)), norm(vecGModelNext) ) );
msg( __FILE__, __LINE__, sprintf( "  ||V^T*g||: (%10.3e) %10.3e -> %10.3e (%10.3e)", norm(matV'*funchGSmooth(vecAIdeal)), norm(matV'*funchGSmooth(vecA)), norm(matV'*funchGSmooth(vecXNext)), norm(matV'*vecGModelNext) ) );

return;
