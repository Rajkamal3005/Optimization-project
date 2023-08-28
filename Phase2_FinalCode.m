

%Group 9/ Raj Kamal Das: 224363005, Vivek Raj: 224363009

clc
file1=fopen("initial_vector.txt",'rt');
x0=textscan(file1,"%f");
fclose(file1);
x0=cell2mat(x0)';






%cauchy's steepest descent
%x0=input("Enter the input vector");
steepest_descent(x0);
function fun=main_objective(x) % multivariable objective function
N=length(x);
f1=0;
for i=1:N
    f1=f1+i*x(i)*x(i);
end
fun=f1;

% -----Question No.2---------
% N=length(x);
% f1=0;
% for i=1:N-1
%     f1=f1+100*(x(i+1)-x(i))^2;
% end
% 
% f2=0;
% for i=1:N-1
%     f2=f2+(x(i)-1)^2;
% end
% fun=f1+f2;
%------Question No.3---------
% N=length(x);
% f1=(x(1)-1)^2;
% f2=0;
% for i=2:N
%     f2=f2+i*(2*x(i)^2-x(i-1))^2;
% end
% fun=f1+f2;
%------Question No.4---------
% N=length(x);
% f1=0;
% for i=1:N
%     f1=f1+(x(i)-1)^2;
% end
% 
% f2=0;
% for i=2:N
%     f2=f2+x(i)*x(i-1);
% end
% fun=f1-f2;
%------Question No.5---------
% N=length(x);
% f1=0;
% for i=1:N
%     f1=f1+x(i)^2;
% end
% 
% f2=0;
% for i=1:N
%     f2=f2+0.5*i*x(i);
% end
% fun=f1+f2^2+f2^4;
end

function steepest_descent(x0) %Function for Steepest descent
%step1
M=100; %Maximum no. of iteration
func_val=[]; % For function value array
f=main_objective(x0); % Fuction value at initial vector
func_val(1)=f;
tol_1=0.0001;  % Tolerance for condition-1
tol_2=0.0001;  % Tolerance for condition-2
k=0;
iteration=[];  % For no. of iterations
iteration(1)=0;


%step-2
grad_val=gradient(x0); %Gradient vector
y=norm(grad_val);
fprintf("The value of y is %.3f :",y);
i=1;
%step-3
while k<M
    if norm(grad_val)<tol_1
        fprintf('The value of x0 is %g\n',x0);
        fprintf("Condition-1");
        fprintf("The minima point lies at %g\n",x0);
        
        break;
    end
    a=-10;b=10; % Lower and upper bound for Exhaustive Search
    val=Exhaustive_Search(a, b,x0,grad_val);
    fprintf("The value of alpha is %8.3f",val);
    xn=x0-val*grad_val; % Updating Initial Vector
    d1=-grad_val/norm(grad_val);
    disp(xn);
    f=main_objective(xn);
    fprintf("The function value is %8.3f\n",f)
    if (norm(xn-x0)/norm(x0))<tol_2
        fprintf("Condition-2\n");
        %fprintf("The minima point lies at %g\n",xn);
        fprintf("...............................\n");
        f=main_objective(xn);
        disp(xn);
        fprintf("..............................\n");
        fprintf("The function value is %8.3f\n",f);
        %condition=false;
        break;
    else
        x0=xn;
        func_val(i+1)=main_objective(x0);
        grad_val=gradient(x0);
        d2=-grad_val/norm(grad_val);
    end
    D=dot(d1,d2);
    Theta = atan2d(1,D); 
    fprintf("The angle is %f\n",Theta);
    if Theta>5 % Condition for linear independance
        fprintf("The directions of vectors are linearly independent.\n");
    end
    k=k+1;
    iteration(i+1)=k;
    fprintf("The value of k is %d\n",k);
    plot(iteration,func_val);
    title("Steepest Descent");
    
    xlabel("No. of iterations");
    ylabel("Function Value");
end
plot(iteration,func_val);
title("Steepest Descent");
    
xlabel("No. of iterations");
ylabel("Function Value");

end

function val=gradient(x0)
N=length(x0);
e=0.001;
grad=zeros(1,N);
for i=1:N
    for j=1:N
        if(i==j)
            x1(j)=x0(j)+e;
        else
            x1(j)=x0(j);
        end
        if(i==j)
            x2(j)=x0(j)-e;
        else
            x2(j)=x0(j);
        end        
        
    end
    f1=main_objective(x1);
    f2=main_objective(x2);
    grad(i)=(f1-f2)/(2*e);
    
end

val=grad;
            
end




%Group 9/ Raj Kamal Das: 224363005, Vivek Raj: 224363009




function fun_val = SingleObjective_Fun(x,x0,grad_val)
    
    xn=x0-x*grad_val;
    fun_val =main_objective(xn);
    
end


function val=Exhaustive_Search(a,b,x0,grad_val)
    fprintf('**EXHAUSTIVE  SEARCH METHOD STARTED**');
    %n = input('\nEnter number of section for exhaustive search (n)  = ');
    n=10;
    % Step 1
    delta = (b - a) / n; % Set data value
    x1 = a; 
    x2 = x1 + delta; 
    x3 = x2 + delta; % New points
    feval = 0; % Number of function evaluations
    f1 = SingleObjective_Fun(x1,x0,grad_val); % Calculate objective function value
    f2 = SingleObjective_Fun(x2,x0,grad_val);
    f3 = SingleObjective_Fun(x3,x0,grad_val);
    fnew=[];
    iteration=[];
    i = 1; % Number of iterations
    feval = feval + 3;
    out = fopen('Project_1.out', 'w'); % Output file
    fprintf(out, 'Outputs for Exhaustive Search Method :\n');
    fprintf(out, '#It\t\t  x1\t\t  x2\t\t  x3\t\t  f(x1)\t\t  f(x2)\t\t  f(x3)\n');
    condition = true;
    while condition
        fprintf(out, '%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n',i,x1,x2,x3,f1,f2,f3);
        % Step 2: Termination condition
        if (f1 >= f2)
            if (f2 <= f3)
                break;
            end
        end
        % If not terminated, update values
        x1 = x2; x2 = x3; x3 = x2 + delta;
        f1 = f2;
        f2 = f3;
        f3 = SingleObjective_Fun(x3,x0,grad_val);
        feval = feval + 1;
        i = i + 1;
        if (x3 > b) % Step 3
            condition = false;
        end
        fnew(length(fnew)+1)=f2;
        iteration(length(iteration)+1)=i;
    end
    
    fprintf('\n*************************\n');
    fprintf('\nThe minimum point lies between (%8.3f, %8.3f)after doing Exhaustive Search Method', x1, x3);
    fprintf('\nTotal number of function evaluations in exhaustive search method : %d\n', feval);
    % Store in the file
    fprintf(out, '\nThe minimum point lies between (%8.3f, %8.3f) after doing Exhaustive Search Method', x1, x3);
    fprintf(out, '\nTotal number of function evaluations in Exhaustive Search Method : %d\n', feval);
    val=bisection(x1,x3,x0,grad_val);

end

function value=bisection(x1,x3,x0,grad_val)
    
    %BISECTION METHOD STARTED

        % x1 and x3 are initial boundary values from Exhaustive Search Method to Bisection Method 
fprintf('\n**BISECTION METHOD STARTED**');
e=0.0001;

i=1;    %No. of iteration for Bisection Method
f1=SingleObjective_Fun(x1+e,x0,grad_val);%To find df(x)/dx
f1=subs(f1);
f2=SingleObjective_Fun(x1-e,x0,grad_val);
f2=subs(f2);
dfx1=(f1-f2)/(2*e);
fprintf("The value is %.3f",dfx1);

f3=SingleObjective_Fun(x3+e,x0,grad_val);
f4=SingleObjective_Fun(x3-e,x0,grad_val);
dfx3=(f3-f4)/(2*e);
feval=4; %function evaluation in Step 1 is 4


%fprintf(out, '\nOutputs for Bisection Method :\n');
%fprintf(out, '#It\t\t  x1\t\t  x3\t\t  c\t\t  f(c)\n');
fnew=[];
iteration=[];
while dfx1*dfx3<0%Step 1       %if any of the value is negative then only executes 
    c=(x1+x3)/2;%Step 2
    fprintf('\n Value of c : %8.3f',c)% middle point c
    fprintf("\tIteration no.= %.d\n",i);
    
    f5=SingleObjective_Fun(c+e,x0,grad_val);
    f5=subs(f5);
    f6=SingleObjective_Fun(c-e,x0,grad_val);
    f6=subs(f6);
    dfc=(f5-f6)/(2*e);
    feval= feval+2; %function evaluation in each iteration is 2 
    
    f7=SingleObjective_Fun(c,x0,grad_val);
    
    %fprintf(out, '%d \t%8.3f \t%8.3f \t%8.3f \t%8.3f  \n',i,x1,x3,c,f7);%for Output file
        
    if abs(dfc)<e       %(Step 3)if absolute value of 1st derivative function value is less than tolerance level then minima is found                      
        fprintf("\nOptimal point lie at x = %.2f\n",c);
        %fprintf(out, '\nOptimal point lie at (%8.3f) after performing Bisection Method\n',c);
        fprintf('\n*****OPTIMIZATION IS OVER*****\n');
        %fprintf(out,'\n\n**OPTIMAL POINT IS REACHED**');
        break
    else          %Updating the boundary values
        if dfc<0
            x1=c;
       
        else
            x3=c;
        end
    end
    i=i+1;
    fnew(length(fnew)+1)=f7;
    iteration(length(iteration)+1)=length(iteration)+i;
    
end

  value=c;  
  fprintf('\n Total number of function evaluation in Bisection Method = %d\n',feval);  
  %fprintf(out,'\nTotal number of function evaluations in Bisection Method : %d\n',feval);  
    
  
  plot(iteration,fnew);
  %axis([0 length(iteration) min(fnew) max(fnew)])
  xlabel('no. of iteration')
  ylabel('function value')
end




