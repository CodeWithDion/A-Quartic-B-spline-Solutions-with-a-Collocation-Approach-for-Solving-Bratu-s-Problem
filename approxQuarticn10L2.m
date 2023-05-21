syms x 
format long
delta_x=0.1;
n=9;
h=0.1;
lamda = 2;
f = @(teta) teta - (sqrt(2.*lamda).*cosh(teta/4));
dfdteta = @(teta) 1 - (sqrt(2.*lamda).*sinh(teta/4))/4;
e = 10.^-4;
teta0 = 0;

%%% Newton Raphson for Calculate the value of Teta in the Exact Solution %%%
for i=1:100
   
    teta = teta0 - f(teta0)/dfdteta(teta0);

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

        S = - 0.204938594900736.*B(1) - 0.060955672758367.*B(2) + 0.064289006091701.*B(3) + 0.168271928234069.*B(4) + 0.248668058261447.*B(5) + 0.303480145433707.*B(6) + 0.331273873192810.*B(7) + 0.331270488608112.*B(8) + 0.303483530018405.*B(9) + 0.248664673676750.*B(10) + 0.168275312818767.*B(11) + 0.064285621507003.*B(12) - 0.060952288173670.*B(13) - 0.204941979485433.*B(14);
  
   
    Error = abs(u-S);
    m=[x(j+2);u;S;Error];
    fprintf('%3.2f   %3.15f   %3.15f   %3.15e\n',m);
end