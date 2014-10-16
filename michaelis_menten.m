function f = michaelis_menten(x,xdata)

f = (x(1).*xdata(:,1))./(x(2)+xdata(:,1));

end

