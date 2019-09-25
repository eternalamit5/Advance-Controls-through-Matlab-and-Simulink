%impresp.m
%calculates the impulse response of a strictly proper, single-input, single-output system

function[y,t]=impresp(num,den,tO,dt,tf);
%Program for calculation of impulse response of strictly proper SISO systems
%num = numerator polynomial coefficients of transfer function
%den = denominator polynomial coefficients of transfer function
%(Coeff icients of 'num' and 'den' are specified as a row vector, in decreasing powers of 's')
%tO = time at which unit impulse input is applied
%dt = time-step (should be smaller than 1/ (largest natural freq.))
%tf = final time for impulse response calculation
%y = impulse response; t= vector of time points
%Find a partial fraction expansion of num/ (den) :-
[r,p,k]=residue(num,den) ;
%Calculate the time points for impulse response : -
t=tO:dt:tf ;
%Find the multiplicity of each pole, p(j):-
for j=1:size(p)
n=1;
    for i=1 :size(p)
        if p(j)==p(i)
            if (i~=j)
            n=n+1 ;
            end
        end
    end
mult(:,j)=n;
end
%Calculate the impulse response by inverse Laplace transform of
%partial-f raction expansion :-
y=zeros(size(t)) ;
j=1;
while j<=size(p,1)
    for i=1:mult(:,j)
    y=y+r(j+i-1)*((t-tO).^(i-1)).*exp(p(j)*(t-tO))/factorial(i-1);
    end
j=j+i;
end