clc;
clear;
close all;
addpath lib

tic;
% -- DATA LOAD ---------------------------------------------------------- %
SC = '27';
LO = '21';
data = load(['sample.mat']).data;%

% -- MAIN --------------------------------------------------------------- %
data_x = data(:,1);
data_y = data(:,2);
data_z = data(:,3);

for n = 1:3
    if n == 1
        NaNx = ~isnan(data_x);
    elseif n == 2
        NaNx = ~isnan(data_y);
    else
        NaNx = ~isnan(data_z);
    end
    data_x = data_x(NaNx);
    data_y = data_y(NaNx);
    data_z = data_z(NaNx);
end
clear data;
toc;disp('DATA LOAD - Done')

% -- RASTERIZATION ------------------------------------------------------ %
MeshSize = 0.05; % [mm]
x_grid = linspace(min(data_x),max(data_x),(max(data_x) - min(data_x))/MeshSize+1);
y_grid = linspace(min(data_y),max(data_y),(max(data_y) - min(data_y))/MeshSize+1);
[X,Y] = meshgrid(x_grid, y_grid);
Z = NaN*ones(size(X));
% zq = griddata(data_x, data_y, data_z, xq, yq);
data = [X(:) Y(:) Z(:)];

[m,~] = size(data);
idx = randperm(m);
data = data(idx,:);

normals = NaN*ones([size(X) 3]);
curvatures = NaN*ones([size(X) 1]);

% -- MOVING LEAST SQUARE & SURFACE NORMAL & CURVATURE ------------------- %
WindowSize = 0.2; % [mm]
[size_x, size_y] = size(X);
for i = 1:size_x
    for j = 1:size_y
        % APPROXIMATION POINT
        xi = [X(i,j), Y(i,j), Z(i,j)];
    
        % POINTS IN SUBTLE AREA
        xyz_idx = ((data_x-xi(1)).^2 + (data_y-xi(2)).^2) <= (WindowSize).^2;
        if length(find(xyz_idx))<2
            X(i,j) = NaN;
            continue
        end
        P_K = [data_x(xyz_idx) data_y(xyz_idx) data_z(xyz_idx)];
    
        % GET Z-VALUES OF THE RECONSTRUCTED SURFACE BY MOVING LEAST SQUARES
        S = P_K;
        [S(:,3),C] = OLS(P_K(:,1),P_K(:,2),P_K(:,3),1);
        f_xi = [xi(1) xi(2) xi(1:2)*C(1:2)+C(3)];
        xi(3) = f_xi(3);
        Z(i,j) = xi(3);
        
        % CALCULATE SURFACE NORMALS
        A = cov(S);
        [V,D] = eig(A); % V are eigenvectors, D are eigenvalues
        [lamda0,idx] = min(diag(D));

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

        % CALCULATE CURVATURES
        curvatures(i,j) = lamda0/sum(diag(D));
        toc;
    end
end
% curvatures = normalize(curvatures);
Z(Z(:) >60) = NaN;
Z(Z(:) <45) = NaN;
toc;disp('RASTERIZATION - Done');

Ra = mean(abs(data_z - OLS(data_x,data_y,data_z,1)));
% 
% %%
% % % -- GABOR FILTER ------------------------------------------------------- %
% % tic;
% % lambda = 10;
% % 
% % B = 1; % bandwidth 
% % gamma = 1; % shape factor of Gaussian fun`ction
% % psi = 0; % offsets
% % 
% % Agray_xy = normals(:,:,1) + 1i*normals(:,:,2);
% % Agray_xy(isnan(Agray_xy)) = 0 + 1i*0;
% % [numRows, numCols] = size(Agray_xy);
% % 
% % deltaTheta = 30;
% % orientation = 0:deltaTheta:(180-deltaTheta);
% % 
% % gb_r = cell(length(lambda)*length(orientation),1);
% % gb_i = cell(length(lambda)*length(orientation),1);
% % gb = cell(length(lambda)*length(orientation),1);
% % gb_features = zeros(numRows, numCols,length(lambda)*length(orientation));
% % 
% % for gb_ld = 1:length(lambda)
% %     for gb_or = 1:length(orientation)
% %         % Filter index
% %         gbi = length(orientation)*(gb_ld-1) + gb_or;
% %         
% %         % calculation for sigma x and y
% %         % -- Sigama calculation based on Bandwidth ---------------------- %
% %         % -- https://doi.org/10.1364/JOSAA.2.001160 --------------------- %
% %         % -- http://www.cs.rug.nl/~imaging/simplecell.html#References --- %
% %         sigma_x = lambda(gb_ld)/pi*sqrt(log(2)/2)*(2^B+1)/(2^B-1);
% %         sigma_y = sigma_x/gamma;
% % 
% %         % define size of the filter
% %         sz = fix(8*max(sigma_y,sigma_x));
% %         if mod(sz,2) == 0 
% %             sz = sz+1;
% %         end
% %         
% %         % creating filter window
% %         [gabor_x, gabor_y]=meshgrid(-fix(sz/2):fix(sz/2),fix(sz/2):-1:fix(-sz/2));
% % 
% %         % Rotation 
% %         x_theta=gabor_x*cos(deg2rad(orientation(gb_or)))+gabor_y*sin(deg2rad(orientation(gb_or)));
% %         y_theta=-gabor_x*sin(deg2rad(orientation(gb_or)))+gabor_y*cos(deg2rad(orientation(gb_or)));
% %             
% %         % gabor filters generation
% %         gb_r{gbi} = exp(-0.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda(gb_ld)*x_theta+psi);
% %         gb_i{gbi} = exp(-0.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*sin(2*pi/lambda(gb_ld)*x_theta+psi);
% %         gb{gbi} = gb_r{gbi} + 1i*gb_i{gbi};
% % 
% %         gb_feature = imfilter(Agray_xy,gb{gbi},'conv');
% %         gb_features(:,:,gbi) = (real(gb_feature).^2 + imag(gb_feature).^2);
% %     end
% % end
% % result_imgs = gb_features;
% % kernels = gb;
% % 
% % reshaped_gbf = reshape(gb_features,numRows*numCols,[]);
% % reshaped_gbf = bsxfun(@minus, reshaped_gbf, mean(reshaped_gbf));
% % reshaped_gbf = bsxfun(@rdivide,reshaped_gbf,std(reshaped_gbf));
% % coeff = pca(reshaped_gbf);
% % feature_xy = reshape(reshaped_gbf*coeff(:,1),numRows,numCols).^2;
% % feature_xy = feature_xy-min(feature_xy,[],'all');
% % feature_xy = feature_xy/max(feature_xy,[],'all');
% % toc;disp('DETECTION PROCESS - Done');
% % 
% % %%% PLOTTING
% % figure(80);imshow(flip(feature_xy,1),[])
% % title('FEATURE XY')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % 
% % figure(81);imshow(imbinarize(flip(feature_xy,1)),[])
% % title('BINARIZED FEATURE XY')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % 
% % %%
% % figure(1)
% % scatter3(data_x,data_y,data_z,5,data_z,'filled');axis tight;axis equal;
% % % plot3(data_x,data_y,data_z,'.');axis equal;axis tight;
% % xlabel('X [mm]')
% % ylabel('Y [mm]')
% % zlabel('Z [mm]')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % view(2)
% % %%
% % figure(2)
% % scatter3(X(:), Y(:), Z(:), 30, Z(:),"filled");axis tight;axis equal;
% % % plot3(X(:), Y(:), Z(:),'.b','MarkerSize',0.5);axis tight;axis equal;
% % xlabel('X [mm]')
% % ylabel('Y [mm]')
% % zlabel('Z [mm]')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % view(2)
% % %%
% % figure(3)
% % N_x = normals(:,:,1);
% % N_y = normals(:,:,2);
% % N_z = normals(:,:,3);
% % N = [N_x(:) N_y(:) N_z(:)];
% % plot3(X(:), Y(:), Z(:),'.b','MarkerSize',.1);hold on;axis tight; axis equal;
% % quiver3(X(:), Y(:), Z(:),N(:,1), N(:,2), N(:,3),'linewidth',1); hold off;
% % % title('SURFACE NORMALS')
% % xlabel('X [mm]')
% % ylabel('Y [mm]')
% % zlabel('Z [mm]')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % % view(90,0)
% % view(2)
% % 
% % %%
% % figure(4)
% % scatter3(X(:), Y(:), curvatures(:), 30, curvatures(:),"filled");axis tight;axis equal;
% % % title('SURFACE NORMALS')
% % xlabel('X [mm]')
% % ylabel('Y [mm]')
% % zlabel('Z [mm]')
% % set(gca,'BoxStyle','full','Box','on','FontSize',32,'FontWeight','bold')
% % view(2)
% % 
