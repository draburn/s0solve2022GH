if (0)
myclear;
setprngstates();
tic
for n=1:9
for t=1:10
	sizeN = 2^n;
	%dat(n,t) = rcond(randn(sizeN,sizeN));
	s = svd(randn(sizeN,sizeN));
	%dat(n,t) = s(1)/s(end-1);
	%dat(n,t) = s(1)/s(end);
	dat(n,t) = s(1);
end
end
toc
end
%
for n=1:7
	%m(n) = median(dat(n,:));
	avg(n) = sum(dat(n,:))/size(dat(:,:),2);
	sqAvg(n) = sum(dat(n,:).^2)/size(dat(:,:),2);
	varSq(n) = sqAvg(n) - (avg(n)^2);
	var(n) = sqrt(varSq(n));
	sizeN(n) = 2^n;
end
plot(sizeN,avg./(sizeN.^0.5),'o-');
grid on;
return;
semilogx( ...
  sizeN, avg+var, '-', ...
  sizeN, avg-var, '-', ...
  sizeN, avg, 'x-' );
grid on;
