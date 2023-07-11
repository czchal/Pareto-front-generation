function [J] = J2(X)
%J2 Summary of this function goes here

g=1;
for i=2:30
    g=g+9/29*X(i);


J= g*(1-sqrt(X(1)/g));

end

