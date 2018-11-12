% Author: Luis Badesa

%%
clear all
close all
clc

%%
disp('Modify:')
disp('- a0: DC term of the Fourier series')
disp('- a(k): expression for sine coeffs')
disp('- b(k): expression for cosine coeffs')
disp('ALL OF THEM IN THE FIRST SECTION OF THE CODE, "Fourier coeffs"')
disp(' ')
disp('FIRST OF ALL:')
disp('Draw the desired frequency respone (first lines of code)')
disp(' ')

%% Desired response
x = -pi:0.01:pi;
y = (x>-pi/2).*(x<pi/2);
figure(1)
plot(x,y)
ylabel('|E|','fontsize',15)
xlabel('\psi','fontsize',15)
grid on
hold on

%% Array characteristics 
prompt = sprintf ('Type: broadside, end-fire array or non? \n(answer "broadside", "end_fire" or "none", as a string in MATLAB) \n');
array_type = input(prompt);
if  ~strcmp(array_type,'broadside') && ~strcmp(array_type,'end_fire') && ~strcmp(array_type,'none')
    error('Wrong type of array introduced.')
end

if  strcmp(array_type,'none')
    prompt = sprintf ('Insert the value of alpha. \n');
    alpha = input(prompt);
    if  ~isscalar(alpha) || (alpha>2*pi) || (alpha<-2*pi)
        error('Insert a value for alpha between -2\pi and 2\pi.')
    end
    
    prompt = sprintf ('What is the value of "d", separation between \neach antenna and the following one in the array? \n(insert a number in terms of lambda, \nEx: introducing 1/2 means d=lambda/2) \n');
    % If given kd, calculate d as d=kd/2*pi.
    d = input(prompt);
    if ~isscalar(d)
        error('Introduce a scalar for the value of "d"')
    end
    
    kd = 2*pi*d;
else
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
end

prompt = sprintf ('What is the number "n" of elements in the array? \n');
n = input(prompt);
if ~isscalar(n) || (n~=floor(n)) || n<=1 || (((n+1)/2)~=floor(((n+1)/2)))
    error('Introduce an odd integer scalar, bigger than 1, for the number of elements in the array (has to be odd because of the Fourier series procedure used)')
end
disp(' ')
disp('NOTE:')
disp('This it the total number of elements in the array, but some of them might not be energized.')

m = (n-1)/2; % m is the number of sine terms in the Fourier series (same as the number of cosine terms)

%% Fourier coeffs

a0=1/2; % Insert the value for a0
for k=1:m
    a(k) = 1/(k*pi)*sin(k*pi/2); % Insert the expression for a_k
    b(k) = 0; % Insert the expression for b_k
end

%% Radiation pattern as a function of psi
psi = -pi:0.001:pi;
sum_coeff =0;
for k=1:m
    sum_coeff = sum_coeff + 2*(a(k)*cos(k*psi)-b(k)*sin(k*psi));
end

E = a0 + sum_coeff;

figure(1)
plot(psi,E)

% % Create the polynomial for the electric field E (see back of page 26 of my notes)
% for k=1:m
%     A(k) = a(m-k+1)-j*b(m-k+1);
% end
% for k=m+2:2*m+1
%     A(k) = a(k-m-1)+j*b(k-m-1);
% end
% A(m+1) = a0;
% 
% psi = -pi:0.001:pi;
% z = exp(j*psi);
% E = 0; % Initialize
% for k=1:m
%     E = E + A(k)*z.^(k-m-1);
% 
% end
% for k=m+2:2*m+1
%     E = E + A(k)*z.^(m-k+1);
% end
% E = E + A(m+1); % This is the z^0 term
% 
% E = abs(E);
%
% E = 1/pi*(-1/3*z.^(-3)+z.^(-1)+pi/2+z-1/3*z.^3);

%% Radiation pattern as a function of phi
clear E

phi = 0:0.001:2*pi;
%psi = (kd*cos(phi)+alpha)
sum_coeff =0;
for k=1:m
    sum_coeff = sum_coeff + 2*(a(k)*cos(k*(kd*cos(phi)+alpha))-b(k)*sin(k*(kd*cos(phi)+alpha)));
end

E = a0 + sum_coeff;
E = abs(E);

figure(2)
polar(phi,E)
