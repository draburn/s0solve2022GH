function datOut = extFit_analyzePt( x0, p, rvecX, rvecF, rvecW=[] )
	thisFile = "extFit_analyzePt";
	msg( thisFile, __LINE__, "This is slower than the original!" );
	%
	% The logic:
	%  d_n = abs( x_n - x_0 )^p
	%  sigma_1  = sum_n w_n
	%  sigma_d  = sum_n w_n * d_n
	%  sigma_dd = sum_n w_n * d_n^2
	%  sigma_f  = sum_n w_n * f_n
	%  sigma_fd = sum_n w_n * f_n * d_n
	%  [ a; b ] = [ sigma_1, sigma_d; sigma_d, sigma_dd ] \ [ sigma_f; sigma_df ]
	%  rho_n = a + b*d_n - f_n
	%  omega = 0.5 * ( sum_n w_n * rho_n^2 ) / sum_1
	%
	if (isempty(rvecW))
		numPts = size(rvecX,2);
		rvecW = ones(1,numPts);
	end
	%
	rvecD = abs(rvecX-x0).^p;
	rvecWD = rvecW.*rvecD;
	rvecWDD = rvecWD.*rvecD;
	rvecWDF = rvecWD.*rvecF;
	rvecWF = rvecW.*rvecF;
	sigma1 = sum(rvecW);
	sigmaF = sum(rvecWF);
	sigmaD  = sum(rvecWD);
	sigmaDD = sum(rvecWDD);
	sigmaDF = sum(rvecWDF);
	%
	denom = sigma1 * sigmaDD - (sigmaD).^2;
	bigA = ( sigmaDD * sigmaF - sigmaD * sigmaDF ) ./ denom;
	bigB = ( sigma1 * sigmaDF - sigmaD * sigmaF ) ./ denom;
	rvecRho = bigA + bigB * rvecD - rvecF;
	rvecWRho = rvecW .* rvecRho;
	rvecWRhoSq = rvecWRho .* rvecRho;
	omega = 0.5 * sum( rvecWRhoSq ) / sigma1;
	%
	datOut.omega = omega;
return;
end

%!test
%!	thisFile = "test extFit_analyzePt_mat 1";
%!	setprngstates();
%!	numFigs = 2;
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
%!	datOut_neo = extFit_analyzePt_mat( y_mat, p_mat, x_rvec, f_rvec );
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
%!	contourf( y_vec, p_vec, datOut_neo.matOmega.^2 );
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
%!	thisFile = "test extFit_analyzePt_mat 1";
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
%!	for i1=1:sizeP
%!	for i2=1:sizeY
%!		datOut_temp = extFit_analyzePt( y_mat(i1,i2), p_mat(i1,i2), x_rvec, f_rvec );
%!		matOmegaNeo(i1,i2) = datOut_temp.omega;
%!	end
%!	end
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
%!	contourf( y_vec, p_vec, matOmegaNeo.^2 );
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