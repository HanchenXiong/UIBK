%*------------------copy all face images into a folder called "all_face" -----------------

mkdir ('all_faces');
dirinfo=dir('orl_faces');

dirinfo(~[dirinfo.isdir])=[];    % remove non-directories


% remove '.' and '..' memebers 
tf=ismember({dirinfo.name}, {'.','..'});    
dirinfo(tf)=[];

numsubdir=size(dirinfo,1);
root_path=pwd;

for i=1:numsubdir
	cd(['orl_faces/',dirinfo(i,1).name]);
	filenames=strcat('*.','pgm');
	filenames=dir(filenames);
	num_files=size(filenames,1);
	for j=1:num_files	
	    im=imread(filenames(j,1).name);
		r_im=imresize(im,0.25);
		imwrite(r_im,[root_path,'/all_faces/s_', num2str(i), '_img_', num2str(j), '.pgm']);
		%copyfile(filenames(j,1).name, [root_path,'/all_faces/s_', num2str(i), '_img_', num2str(j), '.pgm']);
	end;
	cd(root_path);
end;



%*--------------------read all face images and save all them into .mat file--------------
root_path=pwd;
cd('all_faces');
filenames=strcat('*.','pgm');
filenames=dir(filenames);
num_files=size(filenames,1);

im_example=imread(filenames(1,1).name);
im_size=size(im_example);
basis_length=im_size(1)*im_size(2);

data=zeros(basis_length,num_files);

for i=1:num_files
   im_i=imread(filenames(i,1).name);
   r_im_i=reshape(im_i,[basis_length,1]); 
   data(:,i)=r_im_i;
end;

cd(root_path);
save('faces.mat','data');


