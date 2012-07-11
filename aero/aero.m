function [torque] = aero(att_quat, vfs, rho)

if ~exist('rho','var')
    rho = 1.27E-11;
    fprintf('No rho given, assuming 1.27E-11 kg/m^3 (300 km altitude, med-high solar activity)\n');
end

%Define geometry
u = 0.1;        % "1u" = 10 cm
thick=0.005;    % solar panel thickness
%           origin      side x      side y

panels = [  0 0 0       0 u 0       u 0 0 ;     %bottom face
            0 0 3*u     u 0 0       0 u 0 ;     %top face
            0 0 0       0 0 3*u     u 0 0 ;    %XZ face
            0 u 0       0 0 3*u     u 0 0 ;    %other XZ face
            0 0 0       0 0 3*u     0 u 0 ;     %YZ
            u 0 0       0 0 3*u     0 u 0 ;     % other YZ
            0 u 0       u 0 0       0 u 0 ;     %antenna panel
            0 u -thick  u 0 0       0 u 0 ;     %antenna panel
            u 0 0       0 0 3*u     3*u 0 0 ;    %solar panel
            u -thick 0       0 0 3*u     3*u 0 0 ;    %solar panel
            0 0 0       0 0 3*u     -3*u 0 0 ;    %other solar panel
            0 -thick 0       0 0 3*u     -3*u 0 0 ;    %other solar panel
            ];

cg_body = .5 * [ u u 3*u];
        
n_panels=size(panels,1);

%Rotate so freestream is in X
% Tested
vfs_body = quatrotate(att_quat,vfs);
rot_axis = cross(vfs_body/norm(vfs_body),[1 0 0]');
rot_angle = asin(norm(rot_axis));
if(norm(rot_axis)> 1e-10)
    rot_axis = rot_axis/norm(rot_axis);
else
    rot_axis = [0 1 0];
    rot_angle = 0;
end

if (dot(vfs_body,[1 0 0]) < 0)
    rot_angle = rot_angle + pi;
end

q_rot = [cos(rot_angle/2) rot_axis * sin(rot_angle/2)];

cg_fs = quatrotate(q_rot,cg_body);

for i=1:n_panels
    panel_origin(i,:) = quatrotate(q_rot,panels(i,1:3));
    panel_dim1(i,:) = quatrotate(q_rot,panels(i,4:6));
    panel_dim2(i,:) = quatrotate(q_rot,panels(i,7:9));
     
    %find centroid
    panel_centroid(i,:) = panel_origin(i,:) + panel_dim1(i,:)/2 + panel_dim2(i,:)/2;
    
end


q = .5 * rho * norm(vfs)^2;
e_v = [1 0 0];


[discard, index] = sort(panel_centroid(:,1));
index=index';

res = .002;
dim = .5;

foo = @(x)round((x+dim)/res);      %Converts a position to a pixel array index
oof = @(p)res*p-dim;

pixel_panel_index = zeros(2*dim/res, 2*dim/res);
pixel_drag = zeros(2*dim/res, 2*dim/res);
pixel_x = zeros(2*dim/res, 2*dim/res);

Cd=2.5;
drag_per_pixel = .5 * Cd * norm(vfs)^2 * rho * res * res;

for i = 1:length(index)   %iterate over panels, rearmost centroid first
    %e_n = 
    for panel_x = panel_dim1(index(i),:)' * (0:res:1)
        for panel_y = panel_dim2(index(i),:)' * (0:res:1)     %iterate over surface of panel
            panel_pt = panel_origin(index(i),:)' + panel_x + panel_y; %coordinates in fs frame of current panel voxel
            %fprintf('Panel %d, panel_xy = %.01f, %.01f,   
            pixel_pt(1) = foo(panel_pt(2));  %project to freestream YZ plane and turn into an index
            pixel_pt(2) = foo(panel_pt(3));
            pixel_panel_index(pixel_pt(1), pixel_pt(2)) = index(i);     %save the panel number
            pixel_x(pixel_pt(1), pixel_pt(2)) = panel_pt(1);    %save the X depth
            pixel_drag(pixel_pt(1), pixel_pt(2)) = drag_per_pixel;
        end
    end
end

%Torques due to drag
torque_y = 0;
torque_z = 0;
for y=1:size(pixel_drag,1)
    for z=1:size(pixel_drag,2)
        torque_y = torque_y + (oof(y) - cg_fs(2))*pixel_drag(y,z);
        torque_z = torque_z + (oof(z) - cg_fs(3))*pixel_drag(y,z);
    end
end
fprintf('Torques due to drag:\n');
fprintf('\t%e Nm about freestream Y\n',torque_y);
fprintf('\t%e Nm about freestream Z\n',torque_z);


figure(1);
image(flipdim(pixel_panel_index*10,1)); %flip y for plotting
axis square;

%figure(2);
%image(flipdim((pixel_drag>0)*10,1));
%axis square;

fprintf('Total drag = %e N\n',sum(pixel_drag(:)));
%q_rot

