


nx = 10;
ny = 10 ;
ng = nx * ny ;
q = (1:ng)';

a = zeros(ng,1) ;
b = zeros(ng,1) ;
c = zeros(ng,1) ;
d = zeros(ng,1) ;



for i = 1:ng

    point0 = i ; 
    point1 = point0 +1 ;
    point2 = point1 + nx ;
    point3 = point0 + nx ;
    point4 = point3-1;
    point5 = point0-1;
    point6  = point5-nx ;
    point7 = point0 - nx ;
    point8 = point1 - nx ; 

    % 左下角
    if i==1
        d(i) = q(point2);
        % c(i) = q(point4);
        % b(i) = q(point8);
        % a(i) = q(point6);
     
    % 下边
    elseif i>1 && i<nx
        d(i) = q(point2);
        c(i) = q(point4);
        % b(i) = q(point8);
        % a(i) = q(point6);
        
    % 右下角
    elseif i==nx
        % d(i) = q(point2);
        c(i) = q(point4);
        % b(i) = q(point8);
        % a(i) = q(point6);
        
    % 左边
    elseif mod(i,nx)==1 && i>nx && i< ng-nx+1
        d(i) = q(point2);
        % c(i) = q(point4);
        b(i) = q(point8);
        % a(i) = q(point6);
         
    % 右边
    elseif mod(i,nx)==0 && i>nx && i< ng
        % d(i) = q(point2);
        c(i) = q(point4);
        % b(i) = q(point8);
        a(i) = q(point6);
    
    %左上角
elseif i==ng-nx+1
        % d(i) = q(point2);
        % c(i) = q(point4);
        b(i) = q(point8);
        % a(i) = q(point6);
        
    % 上边
    elseif i>ng-nx+1 && i < ng
        % d(i) = q(point2);
        % c(i) = q(point4);
        b(i) = q(point8);
        a(i) = q(point6);
    
    % 右上
    elseif i==ng
        % d(i) = q(point2);
        % c(i) = q(point4);
        % b(i) = q(point8);
        a(i) = q(point6);

    else 
        d(i) = q(point2);
        c(i) = q(point4);
        b(i) = q(point8);
        a(i) = q(point6);
    end

end


ns = ng - (nx+1);
nl = ng - (nx-1);

s =diag(a(end-ns+1:end) , -(nx+1)) + diag(b(end-nl+1:end),-(nx-1)) ...
+ diag(c(1:nl),nx-1) + diag(d(1:ns),nx+1);


spy(s)