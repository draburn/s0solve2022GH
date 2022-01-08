% Sample usage...

function gridImg = xygrids2rgb( gridVX, gridVY, doTranspose=true )
	if (~doTranspose)
		gridVX = gridVX';
		gridVY = gridVY';
	end
	dim1 = size(gridVX,1);
	dim2 = size(gridVX,2);
	assert( isrealarray(gridVX,[dim1,dim2]) );
	assert( isrealarray(gridVY,[dim1,dim2]) );
	scale1 = max([ max(max(abs(gridVX))), max(max(abs(gridVY))) ]);
	scale2 = 5.0;
	scale3 = 10.0;
	gridHSV = double(zeros( dim2, dim1, 3 ));
	gridHSV(:,:,1) = atan2( -gridVY', -gridVX' )/(2.0*pi) + 1.0;gridHSV(:,:,2) = cap( log( 1.0 + (gridVX'.^2 + gridVY'.^2)/10.0 ), 0.5, 1.0 );
	gridHSV(:,:,2) = cap( log( 1.0 + (gridVX'.^2 + gridVY'.^2)/(scale3^2) ), 0.5, 1.0 );
	%gridHSV(:,:,2) = cap( log( 1.0 + (gridVX'.^2 + gridVY'.^2)/10.0 ), 0.5, 1.0 );
	%gridHSV(:,:,3) = cap( sqrt( gridVX'.^2 + gridVY'.^2 ) / scale1, 0.0, 1.0 ).^0.5;
	gridHSV(:,:,3) = cap( asinh(sqrt( gridVX'.^2 + gridVY'.^2 )*scale2) / asinh(scale1/scale2), 0.0, 1.0 );
	gridImg = hsv2rgb(gridHSV);
return;
end

%!test
%!	thisFile = "funcOmega_ellip test: viz";
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	ax = [ -4.0, 4.0, -5.0, 5.0 ];
%!	numXVals = [ 25, 31 ];
%!	x1UniqVals = linspace( ax(1), ax(2), numXVals(1) );
%!	x2UniqVals = linspace( ax(3), ax(4), numXVals(2) );
%!	[ gridX1, gridX2 ] = ndgrid( x1UniqVals, x2UniqVals );
%!	%
%!	gridF = (gridX1).^2 + gridX2.^2;
%!	gridD1 = gridF(3:end,2:end-1) - gridF(1:end-2,2:end-1);
%!	gridD2 = gridF(2:end-1,3:end) - gridF(2:end-1,1:end-2);
%!	gridC1 = gridX1(2:end-1,2:end-1);
%!	gridC2 = gridX2(2:end-1,2:end-1);
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( x1UniqVals, x2UniqVals, gridX1' );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "x1 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( x1UniqVals, x2UniqVals, gridX2' );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "x2 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( x1UniqVals, x2UniqVals, gridF' );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "F vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( x1UniqVals(2:end-1), x2UniqVals(2:end-1), gridD1' );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "dF/dx1 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( x1UniqVals(2:end-1), x2UniqVals(2:end-1), gridD2' );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "dF/dx2 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridC1, gridC2, atan2(gridD2,gridD1) );
%!	colormap(hsv)
%!	axis equal;
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%return
%!	%
%!	numFigs++; figure(numFigs);
%!	image( x1UniqVals(2:end-1), x2UniqVals(2:end-1), xygrids2rgb(gridD1,gridD2) );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "nablaF vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	image( x1UniqVals(2:end-1), x2UniqVals(2:end-1), xygrids2rgb(gridD1*100,gridD2*100) );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "nablaF * 100 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;
%!	%
%!	numFigs++; figure(numFigs);
%!	image( x1UniqVals(2:end-1), x2UniqVals(2:end-1), xygrids2rgb(gridD1/100,gridD2/100) );
%!	axis equal;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "nablaF / 100 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	axis equal;