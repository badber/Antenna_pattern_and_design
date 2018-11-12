clear all
close all
clc

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
if d~=1/2
    disp('The binomial array is not going to be plotted, as "d" has to be 1/2*lambda for this type of array')
end

kd = 2*pi*d;

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

%% Uniform
% Calculate the radiation pattern as a function of phi:
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
        ET_over_E0(i) = n; % Solved by L'Hopitals rule
    else
        ET_over_E0(i) = abs(sin(n*psi/2)/sin(psi/2));
    end
    i=i+1;
end
ET_over_E0 = ET_over_E0/max(ET_over_E0); % Normalize


phi = phi*360/(2*pi); % Convert to degree

figure(1)
plot(phi,ET_over_E0,'Color',[0 0.4470 0.7410])
% Color code:
% 0 0.4470 0.7410
% 0.8500 0.3250 0.0980
% 0.9290 0.6940 0.1250
% 0.4940 0.1840 0.5560
% 0.4660 0.6740 0.1880
% 0.3010 0.7450 0.9330
% 0.6350 0.0780 0.1840
%axis([0 180 0 1])
xlabel('\phi','fontsize',15)
ylabel('E','fontsize',12)
title([num2str(n) '-elements array'])
hold on

clear i psi ET_over_E0

%% Triangular
n = (n+1)/2; % Because the "n" that I use in the expression for the triangular array is not actually its length, it is the length of the uniform array that it comes from (see my class notes back of page 20, the length of the triangular array is 2*n-1, being "n" the length of a uniform array that the triangular array is derived from) 
phi = 0:0.001:pi;
% psi = kd*cos(phi)+alpha;
% z = exp(j*psi);
z = exp(j*(kd*cos(phi)+alpha));
E = 0; % Initialize
for array_position=1:n
    E = E + array_position*z.^(array_position-1);
end
substract = 0; % Initialize. This is used to get the decreasing part of the triangle coefficients
for array_position=(n+1):(2*n-1)
    E = E + (n-1-substract)*z.^(array_position-1);
    substract = substract+1;
end
E = abs(E);
E = E/max(E); % Normalize

phi = phi*360/(2*pi);

figure(1)
plot(phi,E,'Color',[0.8500 0.3250 0.0980])
hold on

%% Binomial
if d==1/2
    % Calculate the basic pattern, the "8" that we get when considering just
    % 2 antennae in the array:
    n_basic = 2;
    
    % Calculate the radiation pattern as a function of phi:
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
    
    ET_over_E0 = ET_over_E0/max(ET_over_E0); % Normalize
    
    phi = phi*360/(2*pi);

    figure(1)
    plot(phi,ET_over_E0,'Color',[0.9290 0.6940 0.1250])
    legend('Uniform','Triangular','Binomial')
end

grid on

