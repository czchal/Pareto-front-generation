function [c,ceq] = con(x,bval,Ubval)
    ceq=[];
    c= [[bval Ubval] + x(31)*n -[J1(x) J2(x)]];

end

