function omega = extFit_analyzePt_quick( x0, p, rvecX, rvecF )
	rvecD = abs(rvecX-x0).^p;
	sigma1 = size(rvecX,2);
	sigmaD = sum(rvecD);
	sigmaDD = sum(rvecD.^2);
	sigmaF = sum(rvecF);
	sigmaDF = sum(rvecD.*rvecF);
	denom = sigma1 * sigmaDD - sigmaD^2;
	bigA = ( sigmaDD * sigmaF  - sigmaD * sigmaDF ) / denom;
	bigB = ( sigma1  * sigmaDF - sigmaD * sigmaF  ) / denom;
	omega = 0.5 * sum( (bigA + bigB*rvecD - rvecF).^2 ) / sigma1;
return;
end

%!test
%!	thisFile = "test extFit_analyzePt_quick 0";
%!	setprngstates(0);
%!	% Pre-load everything...
%!	sizeY = 2;
%!	sizeP = 2;
%!	y_vec = linspace(-3.0,3.0,sizeY);
%!	p_vec = linspace(1.0,5.0,sizeP);
%!	funch_f = @(x)( randn + randn*( abs(x-randn).^(1.0+abs(randn)) ) );
%!	%
%!	[ y_mat, p_mat ] = meshgrid( y_vec, p_vec );
%!	x_rvec = randn(1,5);
%!	f_rvec = funch_f( x_rvec );
%!	%
%!	foo = extFit_analyzePt_quick( y_mat(1,1), p_mat(1,1), x_rvec, f_rvec );
%!	%
%!	foo = 0.5*sum( extFit__getRhoVals( y_mat(1,1), p_mat(1,1), x_rvec, f_rvec ).^2 );
%!	%
%!	foo = extFit_analyzePt( y_mat(1,1), p_mat(1,1), x_rvec, f_rvec );
%!	%
%!	foo = extFit_analyzePt_mat( y_mat, p_mat, x_rvec, f_rvec );


%!test
%!	thisFile = "test extFit_analyzePt_quick 1";
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
%!	x_rvec = randn(1,50000);
%!	f_rvec = funch_f( x_rvec );
%!	%
%!	msg( thisFile, __LINE__, "Doing newnew calc..." );
%!	matOmegaNeo = zeros( sizeP, sizeY );
%!	tic();
%!	for i1=1:sizeP
%!	for i2=1:sizeY
%!		matOmegaNeoNeo(i1,i2) = extFit_analyzePt_quick( y_mat(i1,i2), p_mat(i1,i2), x_rvec, f_rvec );
%!	end
%!	end
%!	toc();
%!	%
%!	msg( thisFile, __LINE__, "Doing new calc..." );
%!	matOmegaNeo = zeros( sizeP, sizeY );
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
%!	msg( thisFile, __LINE__, "Making newnew plot..." );
%!	tic();
%!	numFigs++; figure(numFigs);
%!	contourf( y_vec, p_vec, matOmegaNeoNeo.^2 );
%!	colormap(mycmap(100));
%!	grid on;
%!	toc();
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
