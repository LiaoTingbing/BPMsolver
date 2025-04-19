
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
    m1 = m2 - 1 ;
    m3 = m2 + 1 ;

    % 左边
    if mod(i,nx)==1 
        c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) + q(m3)) + 1/( 0 + q(m2) );
        % a(i) =1/(  q(m1) + q(m2) );
     % 右边
    elseif mod(i,nx)==0    
        % c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) +0) + 1/(  q(m1) + q(m2) );
        a(i) =1/(  q(m1) + q(m2) );
    else
        c(i) = 1/(q(m2) + q(m3));
        b(i) = 1/( q(m2) + q(m3)) + 1/(  q(m1) + q(m2) );
        a(i) =1/(  q(m1) + q(m2) );

    end

end


% ns = ng - (nx+1);
% nl = ng - (nx-1);

s =diag(a(2:end) , -1) + diag(b,0) ...
    + diag(c(1:end-1),1);


spy(s)