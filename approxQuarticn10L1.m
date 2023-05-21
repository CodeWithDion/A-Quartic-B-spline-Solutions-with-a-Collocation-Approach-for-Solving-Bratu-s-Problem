clear all; close all; clc;
syms x 
format longE
delta_x=0.1;
n=9;
h=0.1;
lamda = 1;
f = @(teta) teta - (sqrt(2.*lamda).*cosh(teta/4));
dfdteta = @(teta) 1 - (sqrt(2.*lamda).*sinh(teta/4))/4;
e = 10.^-4;
teta0 = 0;

%%% Newton Raphson for Calculate the value of Teta in the Exact Solution %%%
for i=1:10
   
    teta = teta0 - f(teta0)/dfdteta(teta0)

    if abs(f(teta))<e
        break
    end
    teta0=teta;

end

%%% Calculate the Value of X for every iterations %%%
for i=1:13

    x(i)=(-2+i).*delta_x;

end


%%%% Calculate The Exact and Approximation Solutions %%%%
disp('  x      Exact Solution         S(x)            Error')

for j=1:n

   
    u =-2.*log(cosh((x(j+2)-1/2).*teta/2)/cosh(teta/4));

       B(j) = ((x(j+2)-x(j+2)).^4)/(24*h.^4);
       B(1+j) = ((x(j+2)-x(j+3)).^4)/(24*h.^4);
       B(2+j) = ((x(j+2)-x(j+4)).^4 - 5.*(x(j+2)-x(j+3)).^4)/(24*h.^4);
       B(3+j) = ((x(j+2)-x(j)).^4 - 5.*(x(j+2)-x(j+1)).^4)/(24*h.^4);
       B(4+j) = ((x(j+1)-x(j+2)).^4)/(24*h.^4);
       B(5+j) = ((x(j+2)-x(j+2)).^4)/(24*h.^4);
       B(6+j) = 0;
       B(7+j) = 0;
       B(8+j) = 0;
       B(9+j) = 0;
       B(10+j) = 0;
       B(11+j) = 0;
       B(12+j) = 0;
       B(13+j) = 0;

        S = - 0.091431784616012.*B(1) - 0.026680500175956.*B(2) + 0.028347166842622.*B(3) + 0.073098451282678.*B(4) + 0.107103915456962.*B(5) + 0.129989430272106.*B(6) + 0.141498802356129.*B(7) + 0.141498375780384.*B(8) + 0.129989856847851.*B(9) + 0.107103488881217.*B(10) + 0.073098877858423.*B(11) + 0.028346740266878.*B(12) - 0.026680073600211.*B(13) - 0.091432211191757.*B(14);
  
   
    Error = abs(u-S);
    m=[x(j+2);u;S;Error];
    fprintf('%3.2f   %3.15f   %3.15f   %3.15e\n',m);
end

