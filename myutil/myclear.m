clear;
tic();
%
randState = 0;
%randState = mod(round(time),1E6);
msg("myclear",__LINE__,["randState = " num2str(randState) "."]);
rand("state",randState);
randn("state",randState);
