function datOut = extFit_calcOmega_mat( matBigX, matBigP, rvecX, rvecF, rvecW=[] )
	thisFile = "extFit_calcOmega_mat";
	%
	% The logic:
	%  d_n = abs( x_n - X )^P
	%  sigma_1  = sum_n w_n
	%  sigma_d  = sum_n w_n * d_n
	%  sigma_dd = sum_n w_n * d_n^2
	%  sigma_f  = sum_n w_n * f_n
	%  sigma_fd = sum_n w_n * f_n * d_n
	%  [ A; B ] = [ sigma_1, sigma_d; sigma_d, sigma_dd ] \ [ sigma_f; sigma_df ]
	%  rho_n = A + B*d_n - f_n
	%  omega = 0.5 * ( sum_n w_n * rho_n^2 ) / sum_1
	%
	% But, do this in a way that is quick for a large number of X and P values.
	%
	size1 = size(matBigX,1);
	size2 = size(matBigX,2);
	assert( isrealarray(matBigX,[size1,size2]) );
	assert( isrealarray(matBigP,[size1,size2]) );
	numPts = size(rvecX,2);
	assert( isrealarray(rvecX,[1,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	if (isempty(rvecW))
		rvecW = ones(1,numPts);
	end
	assert( isrealarray(rvecW,[1,numPts]) );
	%
	ary3D   = zeros(size1,size2,numPts);
	ary3WD  = zeros(size1,size2,numPts);
	ary3WDD = zeros(size1,size2,numPts);
	ary3WDF = zeros(size1,size2,numPts);
	parfor n=1:numPts
		ary3D(:,:,n) = abs( rvecX(n) - matBigX ).^matBigP;
		ary3WD(:,:,n) = rvecW(n)*ary3D(:,:,n);
		ary3WDD(:,:,n) = ary3WD(:,:,n).*ary3D(:,:,n);
		ary3WDF(:,:,n) = ary3WD(:,:,n)*rvecF(n);
	end
	rvecWF = rvecW.*rvecF;
	sigma1 = sum(rvecW);
	sigmaF = sum(rvecWF);
	matSigmaD  = sum( ary3WD, 3 );
	matSigmaDD = sum( ary3WDD, 3 );
	matSigmaDF = sum( ary3WDF, 3 );
	%
	matDenom = sigma1 * matSigmaDD - (matSigmaD).^2;
	matBigA = ( matSigmaDD * sigmaF - matSigmaD .* matSigmaDF ) ./ matDenom;
	matBigB = ( sigma1 * matSigmaDF - matSigmaD * sigmaF ) ./ matDenom;
	%
	ary3Rho = zeros(size1,size2,numPts);
	parfor n=1:numPts
		ary3Rho(:,:,n) = matBigA + matBigB .* ary3D(:,:,n) - rvecF(n);
		ary3WRho(:,:,n) = rvecW(n) * ary3Rho(:,:,n);
		ary3WRhoSq(:,:,n) = ary3WRho(:,:,n) .* ary3Rho(:,:,n);
	end
	matOmega = 0.5 * sum( ary3WRhoSq, 3 );
	%
	%
	% Copy main results.
	datOut.matBigA = matBigA;
	datOut.matBigB = matBigB;
	datOut.ary3Rho = ary3Rho;
	datOut.matOmega = matOmega;
	%
	% Copy input.
	datOut.matBigX = matBigX;
	datOut.matBigP = matBigP;
	datOut.rvecX = rvecX;
	datOut.rvecF = rvecF;
	datOut.rvecW = rvecW;
return;
end

%!test
%!	thisFile = "test extFit_calcOmega_mat 1";
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeY = 201;
%!	sizeP = 101;
%!	y_vec = linspace(-3.0,3.0,sizeY);
%!	p_vec = linspace(1.0,5.0,sizeP);
%!	funch_f = @(x)( randn + randn*( abs(x-randn).^(1.0+abs(randn)) ) );
%!	%
%!	[ y_mat, p_mat ] = meshgrid( y_vec, p_vec );
%!	x_rvec = randn(1,5);
%!	f_rvec = funch_f( x_rvec );
%!	%
%!	msg( thisFile, __LINE__, "Doing new calc..." );
%!	tic();
%!	datOut = extFit_calcOmega_mat( y_mat, p_mat, x_rvec, f_rvec );
%!	toc();
%!	%
%!	msg( thisFile, __LINE__, "Doing old calc (Comment out the \"DEPRECATE\" in extFit__getRhoVals.m!)..." );
%!	tic();
%!	omega_mat = zeros(sizeP,sizeY); % Note the "reversed" order of sizes.
%!	for i1=1:sizeP
%!	for i2=1:sizeY
%!		omega_mat(i1,i2) = 0.5*sum( extFit__getRhoVals( y_mat(i1,i2), p_mat(i1,i2), x_rvec, f_rvec ).^2 );
%!	end
%!	end
%!	toc();
%!	%
%!	%
%!	msg( thisFile, __LINE__, "Making new plot..." );
%!	tic();
%!	numFigs++; figure(numFigs);
%!	contourf( y_vec, p_vec, datOut.matOmega.^2 );
%!	colormap(mycmap(100));
%!	grid on;
%!	toc();
%!	%
%!	msg( thisFile, __LINE__, "Making old plot..." );
%!	tic();
%!	numFigs++; figure(numFigs);
%!	contourf( y_vec, p_vec, omega_mat.^2 );
%!	colormap(mycmap(100));
%!	grid on;
%!	toc();



%!test
%!	thisFile = "test extFit_calcOmega_mat 1";
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeY = 2;
%!	sizeP = 2;
%!	y_vec = linspace(-3.0,3.0,sizeY);
%!	p_vec = linspace(1.0,5.0,sizeP);
%!	funch_f = @(x)( randn + randn*( abs(x-randn).^(1.0+abs(randn)) ) );
%!	%
%!	[ y_mat, p_mat ] = meshgrid( y_vec, p_vec );
%!	x_rvec = randn(1,5000);
%!	f_rvec = funch_f( x_rvec );
%!	%
%!	msg( thisFile, __LINE__, "Doing new calc..." );
%!	tic();
%!	datOut = extFit_calcOmega_mat( y_mat, p_mat, x_rvec, f_rvec );
%!	toc();
%!	%
%!	msg( thisFile, __LINE__, "Doing old calc (Comment out the \"DEPRECATE\" in extFit__getRhoVals.m!)..." );
%!	tic();
%!	omega_mat = zeros(sizeP,sizeY); % Note the "reversed" order of sizes.
%!	for i1=1:sizeP
%!	for i2=1:sizeY
%!		omega_mat(i1,i2) = 0.5*sum( extFit__getRhoVals( y_mat(i1,i2), p_mat(i1,i2), x_rvec, f_rvec ).^2 );
%!	end
%!	end
%!	toc();
%!	%
%!	%
%!	msg( thisFile, __LINE__, "Making new plot..." );
%!	tic();
%!	numFigs++; figure(numFigs);
%!	contourf( y_vec, p_vec, datOut.matOmega.^2 );
%!	colormap(mycmap(100));
%!	grid on;
%!	toc();
%!	%
%!	msg( thisFile, __LINE__, "Making old plot..." );
%!	tic();
%!	numFigs++; figure(numFigs);
%!	contourf( y_vec, p_vec, omega_mat.^2 );
%!	colormap(mycmap(100));
%!	grid on;
%!	toc();