tic();
msg("newtoptim20210225_viz",__LINE__,"Generating visualizations...");
numFigs = 0;
%
numFigs++; figure(numFigs);
plot( vecX, vecProbLETau(:,1), 'o-' );
grid on;
hold on;
for n=2:length(list_taus)
	plot( vecX, vecProbLETau(:,n), 'o-' );
end
hold off;
%
numFigs++; figure(numFigs);
plot( vecX, vecF_avgSignPows(:,1), 'o-' );
grid on;
hold on;
for n=2:length(list_pows)
	plot( vecX, vecF_avgSignPows(:,n), 'o-' );
end
hold off;
%
numFigs++; figure(numFigs);
plot( vecX, vecF_avgAbsPows(:,1), 'o-' );
grid on;
hold on;
for n=2:length(list_pows)
	plot( vecX, vecF_avgAbsPows(:,n), 'o-' );
end
hold off;
%
numFigs++; figure(numFigs);
plot( vecX, vecF_signPercentiles(:,1), 'o-' );
grid on;
hold on;
for n=2:length(list_percentiles)
	plot( vecX, vecF_signPercentiles(:,n), 'o-' );
end
hold off;
%
numFigs++; figure(numFigs);
plot( vecX, vecF_absPercentiles(:,1), 'o-' );
grid on;
hold on;
for n=2:length(list_percentiles)
	plot( vecX, vecF_absPercentiles(:,n), 'o-' );
end
hold off;
%
%
toc();
