clear all
close all
clc

prompt = sprintf ('How many arrays are to be compared? \n');
number_arrays = input(prompt);
if ~isscalar(number_arrays) || (number_arrays~=floor(number_arrays)) || number_arrays<=1 || number_arrays>5
    error('Introduce an integer scalar bigger than 1 and smaller than 6 for the number of arrays to be compared')
end

prompt = sprintf ('Broadside or end-fire arrays? \n(answer "broadside" or "end_fire", as a string in MATLAB) \n');
array_type = input(prompt);
if  ~strcmp(array_type,'broadside') && ~strcmp(array_type,'end_fire')
    error('Wrong type of array introduced.')
end

string_array{1} = '1st';
string_array{2} = '2nd';
string_array{3} = '3rd';
string_array{4} = '4th';
string_array{5} = '5th';

for i=1:number_arrays
    prompt = sprintf (['What is the number of elements in the ' string_array{i} ' array? \n']);
    n(i) = input(prompt);
    if ~isscalar(n(i)) || (n(i)~=floor(n(i))) || n(i)<=1
        error('Introduce an integer scalar bigger than 1 for the number of elements in the array')
    end
end

for i=1:number_arrays
    prompt = sprintf (['What is the value of "d", separation between \neach antenna and the following one in the ' string_array{i} ' array? \n(insert a number in terms of lambda, \nEx: introducing 1/2 means d=lambda/2) \n']);
    % If given kd or alpha, do this (see notes of lecture 18, page 4):
    % - If I am given the value of kd, calculate d as d=kd/2*pi.
    % - If I am given the value of alpha and it is different than 0, it means 
    %   that it is an end_fire array, then alpha=kd and d = alpha/2*pi.
    d(i) = input(prompt);
    if ~isscalar(d(i))
        error('Introduce a scalar for the value of "d"')
    end
    if d(i)>=1/2
        disp('NOTE: "d" has to be smaller than 1/2*lambda for nulls to be able to be re-spaced')
    end
    
    kd(i) = 2*pi*d(i);
    
    if  strcmp(array_type,'broadside')
        alpha(i) = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
    elseif strcmp(array_type,'end_fire')
        alpha(i) = -kd(i); % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
    end
end

prompt = sprintf ('Do you want any of these arrays to have uniform current distribution? \n(answer with the number of the array you want to be uniform, \n or 0 if you dont want any array to be uniform) \n');
nulls_respacing = input(prompt);
if  ~isscalar(nulls_respacing) || (nulls_respacing~=floor(nulls_respacing)) || nulls_respacing<0 || nulls_respacing>number_arrays
    error('Insert an integer bigger than 0 and lower than the number of arrays considered')
end
no_respacing(number_arrays) = 0;
if nulls_respacing>0
    no_respacing(nulls_respacing) = 1;
end

%% Patterns
p = 2; % Initialize index
for i=1:number_arrays
    if d(i)>=1/2 || no_respacing(i)==1 % No null re-spacing possible in this case, or not wanted by the user
        % Number of nulls:
        nulls_number(i) = n(i)-1;

        for m=1:n(i)-1
            null{i}{m} = m*(2*pi/n(i));
        end
        null = cell2mat(null{i});
    
        % Draw the radiation pattern as a function of phi:
        phi = 0:0.001:pi;
        %psi = 0:0.0001:pi;
        %psi = kd*cos(phi)+alpha;
        %z = exp(j*psi);
        z = exp(j*(kd(i)*cos(phi)+alpha(i)));
        E = 1; % Initialize
        for q=1:nulls_number(i)
            E = E.*(z-exp(j*null(q)));
        end
        E = abs(E);
        E = E/max(E); % Normalize
        phi = phi*360/(2*pi);

        figure(1)
        plot(phi,E,'LineWidth',2)
        %axis([0 pi 0 n])
        xlabel('\phi','fontsize',18)
        ylabel('|E|','fontsize',15)
        hold on
        
        if d(i)>=1/2
            string_legend{i} = [num2str(n(i)) '-element array, d=\lambda/' num2str(1/d(i)) ' (no null re-spacing possible)'];
        else 
            string_legend{i} = [num2str(n(i)) '-element array, d=\lambda/' num2str(1/d(i)) ', uniform'];
        end
        
        clear null
        
    else % In this case, there is going to be null re-spacing
        %%
        % Calculate the range of psi:
        disp(['Range of psi for the ' string_array{i} ' array:'])
        psi_max(i) = kd(i)+alpha(i); % Using phi=0°
        psi_min(i) = -kd(i)+alpha(i); % Using phi=180°
        psi_max(i)
        psi_min(i)

        % Number of nulls:
        nulls_number(i) = n(i)-1;

        for m=1:n(i)-1
            null{i}{m} = m*(2*pi/n(i));
        end

        figure(p)
        polar(cell2mat(null{i}),ones(1,length(cell2mat(null{i}))),'o')
        hold on
        title([string_array{i} ' array - Original nulls'])

        for l=1:nulls_number(i)
            new_nulls{i}{l} = l*(psi_max(i)-psi_min(i))/nulls_number(i);
        end
        if  strcmp(array_type,'broadside')
            new_nulls = cell2mat(new_nulls{i});
            new_nulls = new_nulls-kd(i)-(psi_max(i)-psi_min(i))/nulls_number(i); 
        elseif strcmp(array_type,'end_fire')
            new_nulls = cell2mat(new_nulls{i});
            new_nulls = new_nulls+(2*pi-2*kd(i))-(psi_max(i)-psi_min(i))/nulls_number(i); 
        end
        p = p+1;
        figure(p)
        polar(new_nulls,ones(1,length(new_nulls)),'o')
        hold on
        title([string_array{i} ' array - New nulls'])
        p = p+1;

        %%
        % Draw the radiation pattern as a function of phi:
        phi = 0:0.001:pi;
        %psi = 0:0.0001:pi;
        %psi = kd*cos(phi)+alpha;
        %z = exp(j*psi);
        z = exp(j*(kd(i)*cos(phi)+alpha(i)));
        E = 1; % Initialize
        for q=1:nulls_number(i)
            E = E.*(z-exp(j*new_nulls(q)));
        end
        E = abs(E);
        E = E/max(E); % Normalize
        phi = phi*360/(2*pi);

        figure(1)
        plot(phi,E,'LineWidth',2)
        %axis([0 pi 0 n])
        xlabel('\phi','fontsize',18)
        ylabel('|E|','fontsize',15)
        hold on
        string_legend{i} = [num2str(n(i)) '-element array, d=\lambda/' num2str(1/d(i)) ', nulls re-spaced'];

        clear new_nulls
    end
end

figure(1)
set(gca,'fontsize',15)
legend(string_legend,'Location','NorthOutside')
grid on
print(figure(1),'-dtiff', 'Radiation_pattern')
