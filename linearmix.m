


function f=linearmix(x,xdata)

    f = x(1) .* xdata(1,:) .* (x(2) + x(3) .* xdata(2,:))+ (1 - x(1)) .* (x(4) + x(5) .* xdata(2,:));

end