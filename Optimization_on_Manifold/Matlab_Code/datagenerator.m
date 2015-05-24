function [pointcloud,T_withou_noise,T] = datagenerator(filename,r1,r2,r3,lx,ly,lz,noise_rate)

%DATAGENERATOR Summary of this function goes here
%   generate expreiment data: point clouds with different motion
pointcloud= load(filename);
num_point=size(pointcloud,1);
pointcloud= pointcloud'; 
%pointcloud= pointcloud*100;    % optional, depends on using metric mm or m 
meanpoint=mean(pointcloud,2);
pointcloud=pointcloud-repmat(meanpoint,1,num_point);  % centralize the point cloud

pointcloud=pointcloud';
save('source.txt', 'pointcloud','-ascii');  

O=pointcloud';
O=[O; ones(1,num_point)];  

% rotation parameters: the angle of yaw, pitch, roll  
phi = r1;     % yaw
chi = r2;     % pitch
psi = r3;     % roll

% tranlation: 
l_x=lx;
l_y=ly;
l_z=lz;

% pose is composed of rotation matrix and translation
pose(1,1)=cos(phi)*cos(chi);
pose(1,2)=cos(phi)*sin(chi)*sin(psi)-sin(phi)*cos(psi);
pose(1,3)=cos(phi)*sin(chi)*cos(psi)+sin(phi)*sin(psi);
pose(1,4)=l_x;
pose(2,1)=sin(phi)*cos(chi);
pose(2,2)=sin(phi)*sin(chi)*sin(psi)+cos(phi)*cos(psi);
pose(2,3)=sin(phi)*sin(chi)*cos(psi)-cos(phi)*sin(psi);
pose(2,4)=l_y;
pose(3,1)=-sin(chi);
pose(3,2)=cos(chi)*sin(psi);
pose(3,3)=cos(chi)*cos(psi);
pose(3,4)=l_z;
pose(4,1)=0;
pose(4,2)=0;
pose(4,3)=0;
pose(4,4)=1;
%disp('rotation:');
%disp(pose(1:3,1:3));
%disp('translation;');
%disp(pose(1:3,4));

% transform original pointcloud with SE(3) pose
T=pose*O;   
T=T(1:3,:);
T=T';
T_withou_noise=T;
save('target.txt','T','-ascii');

noise_mean=zeros(3,num_point);
noise_cov=noise_rate;
noise=normrnd(noise_mean,noise_cov);
T=T+noise';
save('target_noise.txt', 'T','-ascii'); 

end

