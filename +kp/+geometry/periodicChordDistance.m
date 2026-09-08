function [r, delta] = periodicChordDistance(theta, thetaCenters)
%PERIODICCHORDDISTANCE Chord distance on the unit circle.

delta = theta - thetaCenters.';
r = sqrt(2 - 2 * cos(delta));

