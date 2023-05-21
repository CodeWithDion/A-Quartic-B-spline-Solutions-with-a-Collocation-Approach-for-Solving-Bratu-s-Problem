syms x 
format long
delta_x=0.05;
n=19;
h=0.05;
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
for i=1:23

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
       B(14+j) = 0;
       B(15+j) = 0;
       B(16+j) = 0;
       B(17+j) = 0;
       B(18+j) = 0;
       B(19+j) = 0;
       B(20+j) = 0;
       B(21+j) = 0;
       B(22+j) = 0;
       B(23+j) = 0;

       S = - 0.098121532589227.*B(1) - 0.030814882108026.*B(2) + 0.031648215441359.*B(3) + 0.088954865922560.*B(4) + 0.140801093274031.*B(5)...
           + 0.186895617424887.*B(6) + 0.226967266301984.*B(7) + 0.260769105825962.*B(8) + 0.288085774218907.*B(9) + 0.308737050460355.*B(10)...
           + 0.322584119154327.*B(11) + 0.329531574584740.*B(12) + 0.329531679552693.*B(13) + 0.322584014186374.*B(14) + 0.308737155428309.*B(15)...
           + 0.288085669250954.*B(16) + 0.260769210793915.*B(17) + 0.226967161334031.*B(18) + 0.186895722392841.*B(19) + 0.140800988306078.*B(20)...
           + 0.088954970890514.*B(21) + 0.031648110473406.*B(22) - 0.030814777140072.*B(23) - 0.098121637557180.*B(24);
   
    Error = abs(u-S);
    m=[x(j+2);u;S;Error];
    fprintf('%3.2f   %3.15f   %3.15f   %3.15e\n',m);
end
