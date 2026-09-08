function t = wrapPeriodicParameter(t)
%WRAPPERIODICPARAMETER Wrap values to [0,1).

t = mod(t, 1.0);
t(t < 0) = t(t < 0) + 1.0;

