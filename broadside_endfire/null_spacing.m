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
if d>=1/2
    error('"d" has to be smaller than 1/2*lambda for nulls to be able to be re-spaced')
end

kd = 2*pi*d;

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

%%
% Calculate the range of psi:
disp('Range of psi:')
psi_max = kd+alpha % Using phi=0°
psi_min = -kd+alpha % Using phi=180°

% Number of nulls:
nulls_number = n-1;

for m=1:n-1
    null(m) = m*(2*pi/n);
end
figure(1)
polar(null,ones(1,length(null)),'o')
title('Original nulls')

for i=1:nulls_number
    new_nulls(i) = i*(psi_max-psi_min)/nulls_number;
end
if  strcmp(array_type,'broadside')
    new_nulls = new_nulls-kd-(psi_max-psi_min)/nulls_number; 
elseif strcmp(array_type,'end_fire')
    new_nulls = new_nulls+(2*pi-2*kd)-(psi_max-psi_min)/nulls_number; 
end
figure(2)
polar(new_nulls,ones(1,length(null)),'o')
title('New nulls')

%%
% Draw the radiation pattern as a function of phi:
phi = 0:0.001:pi;
%psi = 0:0.0001:pi;
%psi = kd*cos(phi)+alpha;
%z = exp(j*psi);
z = exp(j*(kd*cos(phi)+alpha));
E = 1; % Initialize
for i=1:nulls_number
    E = E.*(z-exp(j*new_nulls(i)));
end
E = abs(E);
E = E/max(E); % Normalize
phi = phi*360/(2*pi);

figure(1)
plot(phi,E)
%axis([0 pi 0 n])
xlabel('\phi','fontsize',15)
ylabel('E_T/E_o','fontsize',12)
title(['n=' num2str(n)])
