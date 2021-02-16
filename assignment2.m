nx = 50;
ny = nx;

G = zeros(nx*ny);
F = zeros(nx*ny,1);
F(1:ny) = 1;
%F((nx-1)*ny+1:nx*ny) = 1;

% Conductivity Map
cMap = ones(ny,nx);
cMap([1:20],[20:30]) = 1e-6;
cMap([30:50],[20:30]) = 1e-6;

ds = 1; % node separation, assumed dx = dy = ds
B = 1 / (ds^2); % coefficient for bulk nodes

for x = 1 : nx
  for y = 1 : ny
    n = (x-1)*ny + y;
    if x == 1
      % left edge
      G(n,n) = 1; % fixed BC
    elseif x == nx
      % right edge
      G(n,n) = 1; % fixed BC
    elseif y == 1
      % bottom edge
      % G(n,n) = 1; % fixed BC
      % free BC (dV/dy = 0)
      G(n, n) = 1;
      G(n, n+1) = -1; 
    elseif y == ny
      % top edge
      % G(n,n) = 1; % fixed BC
      % free BC (dV/dy = 0)
      G(n, n) = 1;
      G(n, n-1) = -1;
    else
      % bulk node
      cxm = (cMap(y,x) + cMap(y,x-1))/2;
      cxp = (cMap(y,x) + cMap(y,x+1))/2;
      cym = (cMap(y,x) + cMap(y-1,x))/2;
      cyp = (cMap(y,x) + cMap(y+1,x))/2;
      G(n,n) = -(cxm + cxp + cym + cyp);
      G(n,n-1) = cym;
      G(n,n+1) = cyp;
      G(n,n-ny) = cxm;
      G(n,n+ny) = cxp;
    end
  end
end

V = G\F;

V_map = reshape(V,[ny,nx]);
[Ex,Ey] = gradient(V_map,ds);
Ex = -1.*Ex;
Ey = -1.*Ey;
Jx = cMap.*Ex;
Jy = cMap.*Ey;
figure(1)
quiver(Jx,Jy);
figure(2)
surf(V_map);