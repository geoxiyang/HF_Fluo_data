
function nee=lightcurve(abc,par)
  
nee = abc(1) -(abc(2)*abc(3).*par)./(abc(2).*par+abc(3));  
