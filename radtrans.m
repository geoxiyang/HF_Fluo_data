


function f = radtrans(x,xdata)
%O2A we use linear for Fmod, linear for rmod
f = (x(1)+xdata(1,:).*x(2)).*xdata(2,:)/pi + (x(3)+xdata(1,:).*x(4));

%O2B we use linear for Fmod, cubic for rmod
% f = (x(1)+xdata(1,:).*x(2)+ xdata(1,:).^2*x(2)+ xdata(1,:).^3*x(3)).*xdata(2,:)/pi + (x(5)+xdata(1,:).*x(6));


end





