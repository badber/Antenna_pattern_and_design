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

kd = 2*pi*d;

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

%%
% Now calculate the radiation pattern as a function of psi:
i=1;
psi = 0:0.001:pi;
ET_over_E0 = [];
for k=1:length(psi)
    if i==1
        ET_over_E0(i) = n; % Solved by L'Hopitals rule
    else
        ET_over_E0(i) = abs(sin(n*psi(k)/2)/sin(psi(k)/2));
    end
    i=i+1;
end
% Now get the rest of the values for psi, since the radiation pattern as a
% function of psi is symmetrical with respect to the psi=pi axis:
psi(i:2*(i-1)) = psi(1:end)+pi;
ET_over_E0(i:2*(i-1)) = fliplr(ET_over_E0);

figure(1)
plot(psi,ET_over_E0)
axis([0 2*pi 0 n])
xlabel('\psi','fontsize',15)
ylabel('E_T/E_o','fontsize',12)
title(['n=' num2str(n)])

clear i psi ET_over_E0

%%
% Now calculate the radiation pattern as a function of phi:
i=1;
phi = [];
ET_over_E0 = [];
if  strcmp(array_type,'broadside')
    % I checked by testing the expression for phi as a function of psi, and this is what I need to print the 360° of phi
    psi_initial = -kd;
    psi_end = kd;
elseif strcmp(array_type,'end_fire')
    psi_initial = -2*kd;
    psi_end = 0;
end
psi_step = (psi_end-psi_initial)/1000;
for psi=psi_initial:psi_step:psi_end
    phi(i) = abs(acos((psi-alpha)/kd));
    if ((psi>-1e-5)&&(psi<1e-5)) || ((psi>2*pi-1e-5)&&(psi<2*pi+1e-5))
        ET_over_E0(i) = n; % Solved by L'Hopitals rule
    else
        ET_over_E0(i) = abs(sin(n*psi/2)/sin(psi/2));
    end
    i=i+1;
end

% Now get the rest of the values for phi, since the radiation pattern as a
% function of phi is symmetrical with respect to the phi=0 axis:
phi(i:2*(i-1)) = -fliplr(phi);
ET_over_E0(i:2*(i-1)) = fliplr(ET_over_E0);

figure(2)
polar(phi,ET_over_E0)
title('E_T/E_o as a function of \phi')

%%
% Calculate width of the primary maxima (Units: radians):
if  strcmp(array_type,'broadside')
    width_PM = 2/(n*d) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
elseif strcmp(array_type,'end_fire')
    width_PM = 2*sqrt(2/(n*d)) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
end

%% THIS PART IS NOT NECESSARY NOW, and actually it doesn't work as the code is now
% % For conversion from psi pattern to phi pattern one also can do this:
% % psi = kd*cos(phi)+alpha, therefore
% % phi = acos((psi-alpha)/kd), but the "acos" function is only defined for
% % values of its argument between -1 and 1, then the values of psi have to
% % be within these limits:
% psi_max = 2*pi*d+alpha
% psi_min = -2*pi*d+alpha
% 
% psi_limits = (psi<=psi_max).*(psi>=psi_min);
% psi_valid = psi.*psi_limits;
% 
% phi = acos((psi-alpha)/kd)
