clear; 
clc;

filename = 'assLegendre.csv';

M = csvread(filename);
[nPts, nMom]  =size(M);

for i =2: nMom
    plot(M(:,1),M(:,i))
    hold on;
end
legend("0,0" , "1,0" , "1,1", "2,0" , "2,1" , "2,2" )