function [r, dotp] = sphereChordDistance(U, V)
%SPHERECHORDDISTANCE Chord distances between points on the unit sphere.

dotp = U * V.';
dotp = min(max(dotp, -1), 1);
r = sqrt(max(2 - 2 * dotp, 0));

