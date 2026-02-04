function [maxvalue, xvalue] = maxqueen(a,b,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


x=linspace(0,5*(a+b),5000);

wvalues=((1-p)*exp(-a./x)+p*exp(-a./x-b./x))./x;

[maxvalue, index]=max(wvalues);
xvalue=x(index);


end

