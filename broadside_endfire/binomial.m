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

d = 1/2; % This means d = lambda/2. For binomial, d must be d<=1/2
kd = pi; % Because kd = 2*pi/lambda*d, and d=lambda/2 for the binomial arrays, so kd = pi

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

% Now calculate the basic pattern, the "8" that we get when considering just
% 2 antennae in the array:
n_basic = 2;

%%
% Now calculate the radiation pattern as a function of psi:
i=1;
psi = 0:0.001:pi;
ET_over_E0 = [];
for k=1:length(psi)
    if i==1
        ET_over_E0(i) = n_basic; % Solved by L'Hopitals rule
    else
        ET_over_E0(i) = abs(sin(n_basic*psi(k)/2)/sin(psi(k)/2));
    end
    i=i+1;
end

% Multiplication by the number of elements in the array, "n":
ET_over_E0 = ET_over_E0.^(n-1);

% Now get the rest of the values for psi, since the radiation pattern as a
% function of psi is symmetrical with respect to the psi=pi axis:
psi(i:2*(i-1)) = psi(1:end)+pi;
ET_over_E0(i:2*(i-1)) = fliplr(ET_over_E0);

psi = psi*360/(2*pi);

figure(1)
plot(psi,ET_over_E0)
axis([0 360 0 n_basic.^(n-1)])
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
    if ((psi>-1e-5)&&(psi<1e-5)) || ((psi>-2*pi-1e-5)&&(psi<-2*pi+1e-5))
        ET_over_E0(i) = n_basic; % Solved by L'Hopitals rule
    else
        ET_over_E0(i) = abs(sin(n_basic*psi/2)/sin(psi/2));
    end
    i=i+1;
end

% Multiplication by the number of elements in the array, "n":
ET_over_E0 = ET_over_E0.^(n-1);

% Now get the rest of the values for phi, since the radiation pattern as a
% function of phi is symmetrical with respect to the phi=0 axis:
phi(i:2*(i-1)) = -fliplr(phi);
ET_over_E0(i:2*(i-1)) = fliplr(ET_over_E0);

figure(2)
polar(phi,ET_over_E0)
title('E_T/E_o as a function of \phi')

% Calculate width of the primary maxima (Units: radian):
if  strcmp(array_type,'broadside')
    width_PM = 2/(n*d) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
    disp('NOTE: this used the formula for uniform distribution, I guess it is also valid for binomial')
elseif strcmp(array_type,'end_fire')
    width_PM = 2*sqrt(2/(n*d)) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
    disp('NOTE: this used the formula for uniform distribution, I guess it is also valid for binomial')
end
