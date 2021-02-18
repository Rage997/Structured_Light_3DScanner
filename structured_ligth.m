% Clear workspace
clc
close all
clear all

%% Given calibration information

K_c = [ 1.4515904453198443e+003, 0., 6.4173428921294862e+002; 0., 1.4543098481380773e+003, 4.6705999995702837e+002; 0., 0., 1. ];

K_p = [ 1.9910689668995669e+003, 0., 5.7522399993777151e+002; 0., 1.9846723239288515e+003, 5.8363251358223056e+002; 0., 0., 1. ];

R_p = [ 9.5017104100048233e-001, 3.6039386628261907e-002, 3.0963874992564855e-001; -3.4715672483391510e-002, 9.9934931692574414e-001, -9.7859513627671715e-003; -3.0978995291662637e-001, -1.4509898370450077e-003, 9.5080391232914441e-001 ];

T_p = [ -2.0069037780485036e+002; -9.1051243942716141e+001;6.0069515527988654e+002 ];

projector_resolution = [1024, 768];

camera_resolution = [1280, 960];

%% load images
ims = cell(18,1);

figure;
% loop through images
for i=1:18
    subplot(4,5,i);
    ims{i} = imread(['big_duck/01/cam_',num2str(i,'%.2d'),'.png']);
    % first image contains color information
    if i==1
        ims{i} = im2double(ims{i});
    end
    % images encoding Gray Code are treated as gray scale
    if i > 2
        ims{i} = rgb2gray(ims{i});
    end
    imshow(ims{i});
end

%% generate Gray Codes

% get size of images that corresponds to camera resolution
[N,M] = size(ims{3});
codes = zeros(N,M);
mask = ones(N,M);
figure;
idx = 1;
% loop through the 8 patterns encoding the Gray Code
for i=3:2:18
    diff = double(ims{i}) - double(ims{i+1});
    % If they do not differ much, discard the point, probably noise
    mask(:,:) = mask & ~(abs(diff) < 255*0.01);
    codes = bitshift(codes, 1);
    codes(:,:) = bitor(codes, diff > 0);
    
    subplot(2,4,idx);
    % display a black and white image for debug
    imshow(diff>0)
    %     imshow(mask)
    idx = idx + 1;
end

%% translate Gray Code to light plane ID
columns = zeros(N,M);

% for each pixel of the camera
for i=1:N
    for j=1:M
        % recover light plane ID from Gray Code
        n=0;
        G = codes(i,j);
        while G > 0
            n = bitxor(n, G);
            G = bitshift(G,-1);
        end
        columns(i,j) = n;
        
    end
end
% mask invalid pixels
columns = mask.*columns;
% Columns relates  camera pixels (i,j) with project columns (plane index)
% (i.e. columns[i,j] = plane_index)
figure;
imagesc(columns); % should smoothly change from blue to red
colormap jet

%% precompute cutting planes

% invert intrinsic projector and camera matrices
T = R_p' * inv(K_p) ;
% precompute plane equations for columns
% store normals and q (parametric equation)
normals = zeros(3, 256); %[n1 n2 n3 ...]
Qs = zeros(3, 256); % [q1 q2 q3 ...]
% consider a 4 bit strip plane (i.e. middle point every 4 bits
for i = 1:projector_resolution(1)/4
    % compute the implicit equation of each light plane
    % defined by the projector location and two rays
    %     P_p = -R_p' * T_p; %projector center
    P_p = -R_p'*T_p;
    ray1 = T*[i*4-2,1, 1]';
    ray2 =  T*[i*4-2, projector_resolution(2), 1]';
    n = cross(ray1, ray2);
    normals(:,i) = n;
    Qs(:,i) = P_p;
end

%% triangulate points
Kc_inv = inv(K_c);

cloud = zeros(1280*960,3)
;
colors = zeros(1280*960,3);

% get pixels of the camera
c = 1:1280;
r = 1:960;
[C,R] = meshgrid(c,r);
C = C(:);
R = R(:);
q = [0;0;0];
ts = zeros(960, 1280);
% loop through all pixels
for i=1:size(C,1)
    % only work on pixels with valid ID to a projected light plane
    plane_idx = columns(R(i),C(i));
    if plane_idx > 0
        % Compute ray from camera r1 = O + lambda d
        d_c = Kc_inv * [C(i), R(i), 1]'; %ray direction from Camera image plane (x,y,1)
        d_c = d_c /norm(d_c);
        % get the implicit light plane equation
        %         for the pixel based on light plane ID
        Q = Qs(:, plane_idx);
        n = normals(:, plane_idx);
        % compute the ray-plane intersection
        lambda = (Q'*n)/ (n'*d_c);
        % check if intersection is valid i.e.
        % intersection point is positive and not at infinity
        max_distance = 600;
        min_distance = 350;
        if lambda > 0 && lambda < max_distance && lambda > min_distance
            % store the new point and add color
            %from first illuminated image
            cloud(i,:)  = lambda*d_c;
            colors(i,:) = ims{1}(R(i),C(i),:);
        end
        
        
    end
end

figure;
scatter3(cloud(:,1), cloud(:,2), cloud(:,3), 10, colors, '.');
axis equal

%% export PLY

ptCloud = pointCloud(cloud, 'Color', colors);
pcwrite(ptCloud,'point_cloud1.ply');