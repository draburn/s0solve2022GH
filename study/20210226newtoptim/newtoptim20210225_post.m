tic();
msg("newtoptim20210225_post",__LINE__,"Performing post-processing...");
%
tauLo = 0.0;
tauHi = 2*absC0 + 3*absC1 + a0 + 3*a1;
if (0)
	list_pows = [ 0.5, 1, 2, 3, 4 ];
	list_percentiles = [ 0.1, 0.3, 0.5, 0.7, 0.9 ];
	list_taus = tauLo + (tauHi-tauLo)*[ 1E-3, 1E-2, 2E-2, 5E-2, 1E-1, 2E-1, 5E-1 ];
else
	list_pows = [2];
	list_percentiles = [ 0.1, 0.5, 0.9 ];
	list_taus = tauLo + (tauHi-tauLo)*[ 1E-3, 1E-2, 1E-1 ];
end
%
%
for n=1:length(list_pows);
	p = list_pows(n);
	%
	matFAbsP = matF_abs.^p;
	vecF_avgAbsPows(:,n) = ( sum( matFAbsP, TRIALS_DIMENSION ) / numTrials ).^(1.0/p);
	%
	matFSignP = matF_sign .* matFAbsP;
	vec_temp = sum( matFSignP, TRIALS_DIMENSION ) / numTrials;
	vecF_avgSignPows(:,n) = sign(vec_temp) .* (abs(vec_temp).^(1.0/p));
	%
	clear vec_temp;
	clear matFAbsP;
	clear matFSignP;
	clear p;
end
for n=1:length(list_percentiles)
	p = list_percentiles(n);
	t = 1 + round((numTrials-1)*p);
	%
	vecF_signPercentiles(:,n) = matF_sort(:,t);
	vecF_absPercentiles(:,n) = matF_absSort(:,t);
	%
	clear t;
	clear p;
end
for n=1:length(list_taus)
	tau = list_taus(n);
	vecProbLETau(:,n) = sum( double(matF_abs<=tau), TRIALS_DIMENSION ) / numTrials;
	clear tau;
end
%
toc();
