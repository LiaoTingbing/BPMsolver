
 
% Exin = h5read("testfile.h5","/Exin") ; 
close;

file = "bpmOut.h5";
x = h5read( file , "/x");
y = h5read( file , "/y");
z = h5read( file , "/z");

index_x = h5read( file , "/index_x");

nx = length(x) ; 
ny = length(y) ;
nz = length(z);


Eout_real = h5read(file,"/Eout_real");
Eout_imag = h5read(file,"/Eout_imag");


Ex  = Eout_real(1:nx*ny,:);
Ex = reshape(Ex , nx,ny,[]);

Ey  = Eout_real(nx*ny+1:end,:);
Ey = reshape(Ey , nx,ny,[]);

for i = 1: nz 
imagesc(Ex(:,:,i))

pause(0.1)
end