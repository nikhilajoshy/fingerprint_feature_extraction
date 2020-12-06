function y=minutie(x) %x is 3*3 window
i=ceil(size(x)/2);    %i=2
if x(i,i)==0; %if x(2,2), ie, center of 3*3 window =0,L=0
    y=0;
else
    y=sum(x(:)) - 1; %else if center =1, L=no:of 1 valued pixels in 3*3 window-1(as center also included)
end