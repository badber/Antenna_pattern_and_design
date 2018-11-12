% Author: Luis Badesa

%%
clear all
close all
clc

%%
prompt = sprintf ('What is the number of elements in the array? \n');
n = input(prompt);
if ~isscalar(n) || (n~=floor(n)) || n<=1
    error('Introduce an integer scalar bigger than 1 for the number of elements in the array')
end

prompt = sprintf ('Broadside or end-fire array? \n(answer "broadside" or "end_fire", as a string in MATLAB) \n');
array_type = input(prompt);
if  ~strcmp(array_type,'broadside') && ~strcmp(array_type,'end_fire')
    error('Wrong type of array introduced.')
end

prompt = sprintf ('What is the value of "d", separation between \neach antenna and the following one in the array? \n(insert a number in terms of lambda, \nEx: introducing 1/2 means d=lambda/2) \n');
% If given kd or alpha, do this (see notes of lecture 18, page 4):
% - If I am given the value of kd, calculate d as d=kd/2*pi.
% - If I am given the value of alpha and it is different than 0, it means 
%   that it is an end_fire array, then alpha=kd and d = alpha/2*pi.
d = input(prompt);
if ~isscalar(d)
    error('Introduce a scalar for the value of "d"')
end

prompt = sprintf ('What is the desired value of the side lobes \nto main lobe attenuation, in dB? \n(introducing "20" means the side lobes is 20dB below the main lobe) \n');
dB_sideLobe = input(prompt);
if ~isscalar(dB_sideLobe) || dB_sideLobe<=0
    error('Introduce a positive scalar')
end

kd = 2*pi*d;

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

%%
Cheby_order = n-1; % The order of the Chebyshev polynomial is one less than the # of elements in the array
% This is valid for up to 40dB, if more attenuation is needed, take more
% points in "x":
if Cheby_order==1
    x = -100:0.01:100;
elseif Cheby_order==2
    x = -10:0.001:10;
else
    x = -3.5:0.001:3.5;
end
y = chebyshevT(Cheby_order,x); 
figure(1)
plot(x,y)
grid on

% Calculate the value of b (see last pages of my notes):
b = 10^(dB_sideLobe/20);

index_zero = find(y>b-1e-2);
x_o = x(index_zero(1));

for k=1:Cheby_order
    delta_k_o(k) = (2*k-1)*pi/(2*Cheby_order);
    x_k_o(k) = cos(delta_k_o(k));
    psi_k_o(k) =2*acos(x_k_o(k)/x_o);
end

phi = 0:0.001:pi;
%psi = 0:0.0001:pi;
%psi = kd*cos(phi)+alpha;
%z = exp(j*psi);
z = exp(j*(kd*cos(phi)+alpha));
E = 1; % Initialize
for l=1:Cheby_order
    E = E.*(z-exp(j*psi_k_o(l)));
end
E = abs(E);
E = E/max(E); % Normalize
phi = phi*360/(2*pi);

figure(1)
plot(phi,E)
%axis([0 pi 0 n])
xlabel('\phi','fontsize',15)
ylabel('|E|','fontsize',12)
grid on
title([num2str(n) '-element array, Tchebysheff distribution'])

