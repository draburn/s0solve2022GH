	clear;
	numFigs = 0;
	setprngstates();
	%setprngstates(10743664); % Image plots glitched. Turn off axis(suqare), "a", zoom out.
	%setprngstates(77077008);
	%setprngstates(5932672)
	%setprngstates(5932672);
	%setprngstates(30660864); % Start on wrong side of a point???
	%setprngstates(58467872);
	%setprngstates(26846592); % Massive slowdown. And, ?!?!?!?!
	%setprngstates(98584480);
	%setprngstates(40121600); % Here, H2 might be better.
	%numPts = 5 + round(abs(randn()*exp(abs(3.0*randn()))))
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
	
	rvecW = ones(size(rvecX));
	prm.numFigs = numFigs;
	extFit_viz( bigX0, bigP0, rvecX, rvecF, [], prm );
