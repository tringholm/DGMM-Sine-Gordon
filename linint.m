function u_new = linint(u,x,x_new)
M = length(u)/2;
u_new = zeros(size(u));
u_new(1) = u(1);
u_new(M) = u(M);
u_new(M+1) = u(M+1);
u_new(end) = u(end);
% for i = 2:M-1
%     k = find(x_new(i) < x,1);
%     x_newdiff = x_new(i) - x(k);
%     x_newdiff2 = x_new(i) - x(k-1)
%     x_olddiff = 1/(x(k) - x(k-1));
%     u_new(i) = (x_newdiff2*u(k) - x_newdiff*u(k-1))*x_olddiff ;
%     u_new(M+i)=(x_newdiff2*u(M+k) - x_newdiff*u(M+k-1))*x_olddiff;
% end
index = zeros(length(2:M-1),1);
for i = 1:length(index)
    index(i) = find(x_new(i+1) < x,1);
end
% index = index(2:end);
x_newdiff = x_new(2:M-1) - x(index);
x_newdiff2 = x_new(2:M-1) - x(index-1);
x_olddiff = 1./(x(index) - x(index-1));
u_new(2:M-1) = (x_newdiff2.*u(index) - x_newdiff.*u(index-1)).*x_olddiff;
u_new(M+2:2*M-1) = (x_newdiff2.*u(M + index) - x_newdiff.*u(M+index-1)).*x_olddiff;


end