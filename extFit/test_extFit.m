	clear;
	thisFile = "test_extFit";
	numFigs = 0;
	%setprngstates();
	%
	setprngstates(33333696); % Stalls
	%setprngstates(55206496); % Slow to converge.
	%setprngstates(40286208); % Yikes; one point is way far away.
	%setprngstates(77173824); % Breaks.
	%
	% DRaburn 2021.06.07: cases below may be boring.
	%setprngstates(19719664); % eye(2,2) may be better here.
	%setprngstates(80489888);
	%setprngstates(10743664); % Image plots glitched. Turn off axis(suqare), "a", zoom out.
	%setprngstates(77077008);
	%setprngstates(5932672)
	%setprngstates(5932672);
	%setprngstates(30660864); % Start on wrong side of a point???
	%setprngstates(58467872);
	%setprngstates(26846592); % Massive slowdown. And, ?!?!?!?!
	%setprngstates(98584480);
	%setprngstates(40121600); % Here, H2 might be better.
	numPts = 5 + round(abs(randn()*exp(abs(randn()))))
	bigX_secret = randn()*exp(abs(3.0*randn()))
	bigP_secret = 1.0 + 3.0*abs(randn())
	bigA_secret = randn()*exp(abs(3.0*randn()));
	bigB_secret = randn()*exp(abs(3.0*randn()));
	rvecX = sort([ ...
	  bigX_secret - abs(randn(1,2)), ...
	  bigX_secret + abs(randn(1,2)), ...
	  bigX_secret + randn(1,numPts-4) ]);
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ).^bigP_secret );
	rvecF = funchF(rvecX);
	rvecW = [];
	prm = [];
	index0 = 1;
	if ( bigB_secret > 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) > rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	elseif ( bigB_secret < 0 )
	while (1)
		if ( (index0==numPts) )
			break;
		elseif ( rvecF(index0+1) < rvecF(index0) )
			break;
		else
			index0++;
			continue;
		end
	end
	end
	bigX0 = (rvecX(index0+1)+rvecX(index0-1))/2.0
	bigP0 = 2.0
	%
	rvecX_first = linspace(min(rvecX),max(rvecX),101);
	rvecF_first = funchF(rvecX_first);
	numFigs++; figure(numFigs);
	plot( ...
	  rvecX, rvecF, 'o', 'markersize', 20, ...
	  rvecX_first, rvecF_first, '*-', ...
	  bigX0*[1,1], [min(rvecF_first),max(rvecF_first)], 'k-' );
	xlabel( "x" );
	ylabel( "f" );
	title( "f vx x" );
	grid on;
	%
	prm_viz.bigX_secret = bigX_secret;
	prm_viz.bigP_secret = bigP_secret;
	prm_viz.numFigs = numFigs;
	extFit_viz( bigX0, bigP0, rvecX, rvecF, [], prm_viz );
	numFigs = 20;
	msg( thisFile, __LINE__, "INITIAL GUESS SHOULD BE UPDATED AFTER EACH POINT! And, should be quad-based." );
	
	bigX = bigX0;
	bigP = bigP0;
	for n=1:10
		[ datFit, retCode ] = extFit( bigX, bigP, rvecX, rvecF );
		%
		bigX = datFit.bigX
		bigP = datFit.bigP
		if ( sum( abs(rvecX-bigX) < eps ) > 0 )
			break;
		end
		rvecX = sort([ rvecX, bigX ]);
		rvecF = funchF( rvecX );
	end
	msg( thisFile, __LINE__, sprintf( "secret values = ( %g, %g ).", bigX_secret, bigP_secret ) );
	msg( thisFile, __LINE__, sprintf( "algori values = ( %g, %g ).", bigX, bigP ) );
	msg( thisFile, __LINE__, sprintf( "residu values = ( %g, %g ).", abs(bigX-bigX_secret), abs(bigP-bigP_secret) ) );
	msg( thisFile, __LINE__, "INITIAL GUESS SHOULD BE UPDATED AFTER EACH POINT! And, should be quad-based." );
	
	prm_viz.bigX_secret = bigX_secret;
	prm_viz.bigP_secret = bigP_secret;
	prm_viz.numFigs = numFigs;
	extFit_viz( bigX, bigP, rvecX, rvecF, [], prm_viz );
	%
	%
	rvecX_fin = linspace(min(rvecX),max(rvecX),101);
	rvecF_fin = funchF(rvecX_fin);
	numFigs++; figure(numFigs);
	plot( ...
	  rvecX, rvecF, 'o', 'markersize', 20, ...
	  rvecX_fin, rvecF_fin, '*-', ...
	  bigX0*[1,1], [min(rvecF_first),max(rvecF_first)], 'y-', ...
	  bigX*[1,1], [min(rvecF_first),max(rvecF_first)], 'k-' );
	xlabel( "x" );
	ylabel( "f" );
	title( "f vx x" );
	axis([0.6,1.1,10.38,10.384])
	grid on;
