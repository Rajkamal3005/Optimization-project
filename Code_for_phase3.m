

%Group 9/ Raj Kamal Das: 224363005, Vivek Raj: 224363009

clc
file1=fopen("initial_vector.txt",'rt');
x0=textscan(file1,"%f");
fclose(file1);
x0=cell2mat(x0)';

% code for phase-3
phase_3(x0);
function [fun,G1,G2]=penalty_fuction(x,R)
g1=(x(1)-5)^2+(x(2)-5)^2-100;
g2=-((x(1)-6)^2+(x(2)-5)^2-82.81);
g3=x(1)-13;
g4=-(x(1)-20);
g5=-(x(2)-4);
f=(x(1)-10)^3+(x(2)-20)^3;
f1=0;

if g1<0
    f1=f1+R*g1^2;
end
if g2<0
    f1=f1+R*g2^2;
end
if g3<0
    f1=f1+R*g3^2;
end
if g4<0
    f1=f1+R*g4^2;
end
if g5<0
    f1=f1+R*g5^2;
end

fun=f1+f;
G1=g1;
G2=g2;
end
function phase_3(x0)
R=0.1;
C=10;
%tol_1=0.001;
tol_2=0.001;
[f,~,~]=penalty_fuction(x0,R);
fprintf("The penalty function value at initial vector is %8.3f",f);
condition=true;
while condition
    x1=steepest_descent(x0,R);
    [fn,g_1,g_2]=penalty_fuction(x1,R);
%     fprintf("The constraint-1 value is %8.3%f",g_1);
%     fprintf("The constraint-1 value is %8.3%f",g_2);
%     if abs(g_1)<tol_1 || abs(g_2)< tol_1
%         fprintf("The minimum vector is\n--------------")
%         disp(xn);
%         fprintf("--------------\n")
%         
%         break;
%     end
    R=C*R;
    xn1=steepest_descent(x1,R);
    [fn1,~,~]=penalty_fuction(xn1,R);
    if abs(fn-fn1)<tol_2
        fprintf("The minimum vector is\n--------------")
        disp(xn1);
        fprintf("--------------\n")
        condition=false;
    end
    x0=xn1;
    
end
end



%cauchy's steepest descent

% function fun=main_objective(x) % multivariable objective function
% N=length(x);
% f1=0;
% for i=1:N
%     f1=f1+i*x(i)*x(i);
% end
% fun=f1;

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
% end

function vector=steepest_descent(x0,R) %Function for Steepest descent
%step1
M=100; %Maximum no. of iteration
func_val=[]; % For function value array
%f=main_objective(x0); % Fuction value at initial vector
[f,~,~]=penalty_fuction(x0,R);
func_val(1)=f;
tol_1=0.0001;  % Tolerance for condition-1
tol_2=0.0001;  % Tolerance for condition-2
k=0;
iteration=[];  % For no. of iterations count
iteration(1)=0;


%step-2
grad_val=gradient(x0,R); %Gradient vector
y=norm(grad_val);
fprintf("The value of y is %.3f :",y);
i=1;
%step-3
while k<M
    if norm(grad_val)<tol_1
        fprintf('The value of x0 is %g\n',x0);
        fprintf("Condition-1");
        fprintf("The minima point lies at %g\n",x0);
        xn=x0;
        break;
    end
    a=-100;b=100; % Lower and upper bound for Exhaustive Search
    val=Exhaustive_Search(a, b,x0,grad_val,R);
    fprintf("The value of alpha is %8.3f",val);
    xn=x0-val*grad_val; % Updating Initial Vector
    d1=-grad_val/norm(grad_val);
    disp(xn);
    %f=main_objective(xn);
    [f,~,~]=penalty_fuction(xn,R);
    fprintf("The function value is %8.3f\n",f)
    if (norm(xn-x0)/norm(x0))<tol_2
        fprintf("Condition-2\n");
        %fprintf("The minima point lies at %g\n",xn);
        fprintf("...............................\n");
        [f,~,~]=penalty_fuction(xn,R);
        disp(xn);
        fprintf("..............................\n");
        fprintf("The function value is %8.3f\n",f);
        %condition=false;
        break;
    else
        x0=xn;
        %func_val(i+1)=main_objective(x0);
        [f1,~,~]=penalty_fuction(x0,R);
        func_val(i+1)=f1;
        grad_val=gradient(x0,R);
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
    i=i+1;
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
vector=xn;

end

function val=gradient(x0,R)
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
    
    [f1,~,~]=penalty_fuction(x1,R);
    [f2,~,~]=penalty_fuction(x2,R);
    grad(i)=(f1-f2)/(2*e);
    
end

val=grad;
            
end




%Group 9/ Raj Kamal Das: 224363005, Vivek Raj: 224363009




function fun_val = SingleObjective_Fun(x,x0,grad_val,R)
    
    xn=x0-x*grad_val;
    %fun_val =main_objective(xn);
    [f1,~,~]=penalty_fuction(xn,R);
    fun_val =f1;
    
end


function val=Exhaustive_Search(a,b,x0,grad_val,R)
    fprintf('**EXHAUSTIVE  SEARCH METHOD STARTED**');
    %n = input('\nEnter number of section for exhaustive search (n)  = ');
    n=100;
    % Step 1
    delta = (b - a) / n; % Set data value
    x1 = a; 
    x2 = x1 + delta; 
    x3 = x2 + delta; % New points
    feval = 0; % Number of function evaluations
    f1 = SingleObjective_Fun(x1,x0,grad_val,R); % Calculate objective function value
    f2 = SingleObjective_Fun(x2,x0,grad_val,R);
    f3 = SingleObjective_Fun(x3,x0,grad_val,R);
    fnew=[];
    iteration=[];
    i = 1; % Number of iterations
    feval = feval + 3;
    %out = fopen('Project_1.out', 'w'); % Output file
    %fprintf(out, 'Outputs for Exhaustive Search Method :\n');
    %fprintf(out, '#It\t\t  x1\t\t  x2\t\t  x3\t\t  f(x1)\t\t  f(x2)\t\t  f(x3)\n');
    condition = true;
    while condition
        %fprintf(out, '%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n',i,x1,x2,x3,f1,f2,f3);
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
        f3 = SingleObjective_Fun(x3,x0,grad_val,R);
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
    %fprintf(out, '\nThe minimum point lies between (%8.3f, %8.3f) after doing Exhaustive Search Method', x1, x3);
    %fprintf(out, '\nTotal number of function evaluations in Exhaustive Search Method : %d\n', feval);
    val=bisection(x1,x3,x0,grad_val,R);

end

function value=bisection(x1,x3,x0,grad_val,R)
    
    %BISECTION METHOD STARTED

        % x1 and x3 are initial boundary values from Exhaustive Search Method to Bisection Method 
fprintf('\n**BISECTION METHOD STARTED**');
e=0.001;

i=1;    %No. of iteration for Bisection Method
f1=SingleObjective_Fun(x1+e,x0,grad_val,R);%To find df(x)/dx

f2=SingleObjective_Fun(x1-e,x0,grad_val,R);

dfx1=(f1-f2)/(2*e);
fprintf("The value of f1 is %.3f",f1);
fprintf("The value of f2 is %.3f",f2);
fprintf("The value of dfx1 is %.3f",dfx1);

f3=SingleObjective_Fun(x3+e,x0,grad_val,R);
f4=SingleObjective_Fun(x3-e,x0,grad_val,R);
dfx3=(f3-f4)/(2*e);
feval=4; %function evaluation in Step 1 is 4
fprintf("The value of dfx3 is %.3f",dfx3);

%fprintf(out, '\nOutputs for Bisection Method :\n');
%fprintf(out, '#It\t\t  x1\t\t  x3\t\t  c\t\t  f(c)\n');
fnew=[];
iteration=[];
% if abs(x3-x1)<e
%     c=x1;
%    
% end
%c=0;
%c=(x1+x3)/2;
while dfx1*dfx3<0 %&& dfx3>0 && dfx1<0%Step 1       %if any of the value is negative then only executes 
    c=(x1+x3)/2;%Step 2
    fprintf("\n value of x1 %d",x1);
    fprintf("\n value of x3 %d",x3);
    
    fprintf('\n Value of c : %8.3f',c);% middle point c
    fprintf("\tIteration no.= %.d\n",i);
    %%fprintf("value of dfx1 %8.3f",dfx1);
    %%fprintf("value of dfx3 %8.3f",dfx3);
    
    f5=SingleObjective_Fun(c+e,x0,grad_val,R);
    %%fprintf(" okkkkkk");
    f6=SingleObjective_Fun(c-e,x0,grad_val,R);
    
    dfc=(f5-f6)/(2*e);
    feval= feval+2; %function evaluation in each iteration is 2 
    fprintf("value of func %8.3f",dfc);
    f7=SingleObjective_Fun(c,x0,grad_val,R);
    
    %fprintf(out, '%d \t%8.3f \t%8.3f \t%8.3f \t%8.3f  \n',i,x1,x3,c,f7);%for Output file
    fprintf("value of dfc %8.3f",dfc);
    if abs(dfc)<e       %(Step 3)if absolute value of 1st derivative function value is less than tolerance level then minima is found                      
        fprintf("\nOptimal point lie at x = %.2f\n",c);
        %fprintf(out, '\nOptimal point lie at (%8.3f) after performing Bisection Method\n',c);
        fprintf('\n*****OPTIMIZATION IS OVER*****\n');
        %fprintf(out,'\n\n**OPTIMAL POINT IS REACHED**');
        break
    else          %Updating the boundary values
        if dfc<0
            x1=c;
            %fprintf("ok222222");
       
        else
            x3=c;
            fprintf("value of x3 %d \n",x3);
            fprintf("value of x1 %d \n",x1);
            fprintf("ok33333\n");
        end

    end

    i=i+1;
    fnew(length(fnew)+1)=f7;
    iteration(length(iteration)+1)=length(iteration)+i;
    
end
  
    
  fprintf('\n Total number of function evaluation in Bisection Method = %d\n',feval);  
  %fprintf(out,'\nTotal number of function evaluations in Bisection Method : %d\n',feval);  
  value=c;  
  
  plot(iteration,fnew);
  %axis([0 length(iteration) min(fnew) max(fnew)])
  xlabel('no. of iteration')
  ylabel('function value')
end




