
clear;

nx = 10;
ny = 10 ;
ng = nx * ny ;
q = (1:ng)';

a = zeros(ng,1) ;
b = zeros(ng,1) ;
c = zeros(ng,1) ;
d = zeros(ng,1) ;

for i = 1:ng
    m2 = i;
    m1 = m2 - nx ;
    m3 = m2 + nx ;
    % 下
    if i>=1 && i<=nx 
        c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) + q(m3)) + 1/( 0 + q(m2) );
        % a(i) =1/(  q(m1) + q(m2) );
     % 上边
    elseif i>=ng-nx+1 && i<=ng    
        % c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) +0) + 1/(  q(m1) + q(m2) );
        a(i) =1/(  q(m1) + q(m2) );
    else
        c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) + q(m3)) + 1/(  q(m1) + q(m2) );
        a(i) =1/(  q(m1) + q(m2) );
    end
end
s =diag(a(nx+1:end) , -nx) + diag(b,0) ...
    + diag(c(1:end-nx),nx);
spy(s)