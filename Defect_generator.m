clc;
clear;
close all;

%%%%%%%%%%%%% INSTRUCTION %%%%%%%%%%%%%%%%%
% SELECT YOUR DEFECT TYPE BELOW and RUN: 
% YOUR OPTIONS -> BUMP, DENT, SCRATCH, BURR
% THEN RUN "RECONSTRUCTION" line 77
% THEN RUN "SURFACE NORMAL" line 126
%% -- BUMP ------------------------------------------------------------- %%
szgrid = 1000;
[xx, yy] = meshgrid(linspace(-5,5,szgrid),linspace(-5,5,szgrid));
zz = exp(-xx.^2-yy.^2) +1;
xx = xx + 5;
yy = yy + 5;

%% OLS 函數定義
function [output, C] = OLS(x, y, z, flag)
    A = [x, y, ones(size(x))];
    C = A \ z;
    output = [x, y, A * C];
end

% %% -- DENT ------------------------------------------------------------- %%
% szgrid = 1000;
% [xx, yy] = meshgrid(linspace(-5,5,szgrid),linspace(-5,5,szgrid));
% zz = -exp(-xx.^2-yy.^2) +2;
% xx = xx + 5;
% yy = yy + 5;
% %% -- SCRATCH ---------------------------------------------------------- %%
% szgrid = 5000;
% [xx, yy] = meshgrid(linspace(-5,5,szgrid),linspace(-5,5,szgrid));
% zz = -exp(-xx.^2-yy.^2) +2;
% z = zz;
% % for i = 1:size(zz)
% %     rnd = randi([-5,5]);
% %     zz(i,6+rnd:95+rnd) = zz(51,6:95);
% % end
% for i = 1:size(zz,1)
%     rnd = randi([-20,20]);
%     if i <= 70
%         zz(i,21+rnd:end-20+rnd) = z(szgrid/2-70+i,21:end-20);
%     elseif i > szgrid-70
%         zz(i,21+rnd:end-20+rnd) = z(i-180,21:end-20);
%         disp(i)
%     else
%         zz(i,21+rnd:end-20+rnd) = z(szgrid/2+1,21:end-20);
%     end
% end
% rnd = randi([0,90]);
% xx_prime = xx*cos(rnd)+yy*sin(rnd);
% yy_prime = -xx*sin(rnd)+yy*cos(rnd);
% xx = xx_prime;
% yy = yy_prime;
% xx = xx + 5;
% yy = yy + 5;
% %% -- BURR ------------------------------------------------------------- %%
% szgrid = 1000;
% [xx, yy] = meshgrid(linspace(-5,5,szgrid),linspace(-0,10,szgrid));
% zz_BUMP = exp(-xx.^2-yy.^2) +1;
% zz_DENT = -exp(-xx.^2-yy.^2) +1;
% % zz =[zz_BUMP(:,26:75),zz_DENT(:,26:75)];
% zz = [zz_BUMP(:,26:75),zz_DENT(:,76:100),zz_DENT(:,76:100)];
% xx = xx + 5;
% surf(xx,yy,zz);
% % shading interp;
% % view(-36,56)
%% RECONSTRUCTION
MeshSize = 0.05; % [mm]
x_grid = linspace(min(min(xx)),max(max(xx)),(max(max(xx))-min(min(xx)))/MeshSize+1);
y_grid = linspace(min(min(yy)),max(max(yy)),(max(max(yy))-min(min(yy)))/MeshSize+1);
[X,Y] = meshgrid(x_grid, y_grid);
Z = NaN*ones(size(X));
normals = NaN*ones([size(X) 3]);

WindowSize = 1*MeshSize; % [mm]
[size_x, size_y] = size(X);
for i = 1:size_x
    for j = 1:size_y
        % APPROXIMATION POINT
        xi = [X(i,j), Y(i,j), Z(i,j)];
    
        % POINTS IN SUBTLE AREA
        xyz_idx = ((xx-xi(1)).^2 + (yy-xi(2)).^2) <= (WindowSize).^2;
        if length(find(xyz_idx))<2
            X(i,j) = NaN;
            continue
        end
        P_K = [xx(xyz_idx) yy(xyz_idx) zz(xyz_idx)];
    
        % GET Z-VALUES OF THE RECONSTRUCTED SURFACE BY MOVING LEAST SQUARES
        S = P_K;
%       [S(:,3),C] = OLS(P_K(:,1),P_K(:,2),P_K(:,3),1);
        [S(:,3), C] = OLS(P_K(:,1), P_K(:,2), P_K(:,3), 1);
        f_xi = [xi(1) xi(2) xi(1:2)*C(1:2)+C(3)];
        xi(3) = f_xi(3);
        Z(i,j) = xi(3) + sum(sum((f_xi - xi).^2)*WeightFunction(S,xi,MeshSize));
%         Z(i,j) = xi(3) + mean(mean((f_xi - xi).^2)*WeightFunction(S,xi,MeshSize));
        % CALCULATE SURFACE NORMALS
%         S = [S; xi];
        A = cov(S);
        [V,D] = eig(A); % V are eigenvectors, D are eigenvalues
        [~,idx] = min(diag(D));

        temp_normal = V(:,idx)./norm(V(:,idx));
        SensorCenter = [X(i,j), Y(i,j), 60];
        p1 = SensorCenter - xi;
        p2 = temp_normal';
        % Flip the normal vector if it is not pointing towards the sensor.
        angle = atan2(norm(cross(p1,p2)),p1*p2');
        if angle > pi/2 || angle < -pi/2
            temp_normal = -temp_normal;
        end
        normals(i,j,:) = temp_normal;
    end
end



%% SURFACE NORMAL (TOP)
N_x = normals(:,:,1);
N_y = normals(:,:,2);
N_z = normals(:,:,3);
% N_z(N_z <0) = 0.7;

figure(3)
N = [N_x(:) N_y(:) N_z(:)];
plot3(X(:), Y(:), Z(:),'.b');hold on;axis tight; axis equal;
% plot3(P_K(:,1),P_K(:,2),P_K(:,3),'.r')
% plot3(S(:,1),S(:,2),S(:,3),'.k')
% plot3(xi(:,1),xi(:,2),xi(:,3),'.k','MarkerSize',30)
quiver3(X(:), Y(:), Z(:),N(:,1), N(:,2), N(:,3)); hold off;
set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
view(2)
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

figure(4);
surf(xx,yy,zz);axis equal;
set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
view(2)
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
shading interp;

data = [X(:), Y(:), Z(:)];
plot3(data(:,1),data(:,2),data(:,3),'.')
%% SURFACE NORMAL (SIDE)
N_x = normals(:,:,1);
N_y = normals(:,:,2);
N_z = normals(:,:,3);
% N_z(N_z <0) = 0.7;

figure(3)
N = [N_x(:) N_y(:) N_z(:)];
plot3(X(:), Y(:), Z(:),'.b');hold on;axis tight; axis equal;
% plot3(P_K(:,1),P_K(:,2),P_K(:,3),'.r')
% plot3(S(:,1),S(:,2),S(:,3),'.k')
% plot3(xi(:,1),xi(:,2),xi(:,3),'.k','MarkerSize',30)
quiver3(X(:), Y(:), Z(:),N(:,1), N(:,2), N(:,3)); hold off;
set(gca,'BoxStyle','full','Box','on','FontSize',64,'FontWeight','bold')
view(0,0)
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

figure(4);
surf(xx,yy,zz);axis equal;
set(gca,'BoxStyle','full','Box','on','FontSize',64,'FontWeight','bold')
view(0,0)
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
shading interp;

