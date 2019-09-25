%Function stepresp.m, 
%calculates the step response of a proper, single-input, single-output system

function [y,t]=stepresp(num,den,t0,dt,tf);
%Program for calculation of step response of proper SISO systems
%num = numerator polynomial coefficients of transfer function
%den = denominator polynomial coefficients of transfer function
%(Coefficients of 'num' and 'den' are specified as a row vector, in
%decreasing powers of 's')
%tO = time at which unit step input is applied
%dt = time-step (should be smaller than 1/(largest natural freq.))
%tf = final time for step response calculation
%y = step response; t= vector of time points
%Find a partial fraction expansion of num/(den.s):-
[r,p,k]=residue(num,conv(den,[1 0]));
%Calculate the time points for step response:-
t=t0:dt:tf;
%Find the multiplicity of each pole, p ( j ) : -
for j=1:size(p)
n=1;
    for i=1:size(p)
        if p(j)==p(i)
            if(i~=j)
                n=n+1;
            end
        end
    end
mult(:,j)=n;
end
%Calculate the step response by inverse Laplace transform of
%partial-fraction expansion:-
y=zeros(size(t));
j=1;
while j<=size(p,1)
    for i=1:mult(:,j)
    y=y+r(j+i-1)*((t-t0).^(i-1)).*exp(p(j)*(t-t0))/factorial(i-l);
    end
    j=j+i;
end