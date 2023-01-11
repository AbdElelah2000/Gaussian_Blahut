%Given peak power constraints of 0.1, 1 and 10
A = [0.1, 1, 10];

%Repeat for each A value
for i = A
    range_of_x = -i:(2*i/1000):i;% Consider a sampling rate of 1000 samples
    p_y_x = [];% Start with empty prob of y given x
    for mean = range_of_x %compute p_y_x
       row_of_y = normpdf(range_of_x,mean,1); % Discretize the output prob mat by setting the mean as each of the samples
       p_y_x( end+1, : ) = row_of_y; % Add it to the p_y_x mat
    end
    % Call the Blahut Arimoto Algorithm function and plot it for each A value
    [Capacity_Val, p_x_star] = blahut_Algo_fn(p_y_x, i);
    Capacity_Val %Print the capacity value in the console
    %Graph the capacity acheiving dist for each A value
    figure
    plot(range_of_x,p_x_star ,'color','red')
    title(sprintf('The capacity-achieving distribution for A = %1.1f', i))
    ylabel('p*(x)')
    xlabel('x')
    legend('AbdElelah Arafah, 400197623')
end


%define a function to modularize each component of the code
function [Capac, p_x] = blahut_Algo_fn(p_y_x_arg, A)
fprintf('For A = %1.1f', A);

%throw an error if theres a column of 0 in p_y_x and set the capacity to 0
column_sum = sum(p_y_x_arg);
if ~isempty(find(column_sum == 0, 1))
    disp('Error: 0 column in p_y_x');
    Capac = 0.0; return;
end

% Size of the sample, and the output
[size_of_x, size_of_y] = size(p_y_x_arg);

%Define and discretize the initial uniform input distribution
p_x = ones(1,size_of_x)/size_of_x;
%initialize q_x with 0s
q_x = zeros(size_of_x,size_of_y);
%Allow an error of 1e-12 for the stopping criterion as recommended by online resources
threshold_err = 1e-12/size_of_x;
%The next increment of p_x in the algo
p_x_increment = [];

%maximizing p_x and q_x for 7500 iterations
for count = 1:size_of_x
    p_y_x_arg(count,:) = p_y_x_arg(count,:)/sum(p_y_x_arg(count,:));
end
for num_of_iterations = 1:7500
    %compute q_x
    for j = 1:size_of_y
        q_x(:,j) = p_x'.*p_y_x_arg(:,j);
        q_x(:,j) = q_x(:,j)/sum(q_x(:,j));
    end
    
    %compute p_x
    for count = 1:size_of_x
        p_x_increment(count) = p_x(count)*prod((p_y_x_arg(count,:)./sum(p_x'.*p_y_x_arg(:,count))).^p_y_x_arg(count,:));
    end
    p_x_increment = p_x_increment/sum(p_x_increment);

    %check if p_x is optimal
    if norm(p_x_increment - p_x) < threshold_err
        break
    else
        p_x = p_x_increment;
    end
end

%initialize the capacity with 0
Capac = 0;

%Compute the capacity for p_x values that are positive
for count = 1:size_of_x
    for j = 1:size_of_y
        if p_x(count) > 0 && q_x(count,j) > 0 %For positive values of q_x and p_x add to the capacity
            Capac = Capac + p_x(count)*p_y_x_arg(count,j)* log2(q_x(count,j)/p_x(count)); %Compute the capacity with log base 2
        end
    end
end
end


