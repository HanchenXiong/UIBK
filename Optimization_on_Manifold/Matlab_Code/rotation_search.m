function [optimal_rotation, trace]=rotation_search(source,target,method)


%---------------when the number points is too large, only subsamples are used---------------
num_point=size(source,1);
rand_idx=randperm(num_point);

sub_sample_perc=0.2;
num_sample=round(sub_sample_perc*num_point);
rand_idx=rand_idx(1:num_sample);       % only sub samples are used to compute objectives 

s_source=source(rand_idx,:);


%-------------------------------basic set up -----------------------------------------------
n=1000;      % number of particles  or SA
T=1540;      % target temperature
delta=0.15;  % opimal delta
nu=2.38/3;
gamma=0.2;
power=1.3;
ESS_threshold=0.5;
size_subsample=1000;

trace=[];    % record the best resutls out of the particle population

%-------------------------------------------------------------------------------------------
if strcmp(method,'N-SA')	
	%--------------initilization of N SA--------------------
%	initials=zeros(9,n);
%	for i=1:n
%		phi = rand*2*pi;   % yaw
%        chi = rand*pi;     % pitch
%        psi = rand*pi;     % roll
%        
%		pose=zeros(3,3);
%		pose(1,1)=cos(phi)*cos(chi);
%		pose(1,2)=cos(phi)*sin(chi)*sin(psi)-sin(phi)*cos(psi);
%		pose(1,3)=cos(phi)*sin(chi)*cos(psi)+sin(phi)*sin(psi);
%		pose(2,1)=sin(phi)*cos(chi);
%		pose(2,2)=sin(phi)*sin(chi)*sin(psi)+cos(phi)*cos(psi);
%		pose(2,3)=sin(phi)*sin(chi)*cos(psi)-cos(phi)*sin(psi);
%		pose(3,1)=-sin(chi);
%		pose(3,2)=cos(chi)*sin(psi);
%		pose(3,3)=cos(chi)*cos(psi);
% 
%		initials(:,i)=reshape(pose,9,1);
%	end;
	initials=rand(9,n)*2-1;
    for i=1:n
		i_initial=initials(:,i);
		i_matrix=reshape(i_initial,3,3);
	    [U,sigma,V]=svd(i_matrix);
		i_rotation=U*V';
		initials(:,i)=reshape(i_rotation,9,1);
	end;
	current_states=initials;
       	
	%-------------------------------------------------------
	k=0;	
	t=1;
	while t<T
		k=k+1
		t=k^power;   % temperature at each iteration
		global_pro=1/(1+k^0.1);
       
        max_objective=-inf; 
	    for i=1:n	 	
		    i_current_state=current_states(:,i);	
			i_current_rotation=reshape(i_current_state,3,3);        
			i_current_objective=computing_objective(s_source,target,i_current_rotation);	
			if i_current_objective>max_objective
				max_objective=i_current_objective;
				optimal_rotation=i_current_rotation;
			end;

			r=rand;
			if (r<global_pro)
				% go to global exploration
				i_next_state=rand(9,1)*2-1;
			else
				% go to global exploration
				i_next_state=mvnrnd(i_current_state, (nu^2)*eye(9));
			end;

			i_next_matrix=reshape(i_next_state,3,3);
			[U,sigma,V]=svd(i_next_matrix);
			i_next_rotation=U*V';
			i_next_objective=computing_objective(s_source, target, i_next_rotation);

	        prob=((i_next_objective+delta)/(i_current_objective+delta))^t;
			
			r=rand;
			if r<prob 
	%			disp('accept');
				i_next_state=reshape(i_next_rotation,9,1);
				current_states(:,i)=i_next_state;
			end;
		end; 

		max_objective
		trace=[trace,max_objective];
	end;    
elseif (strcmp(method,'SMC'))
    %--------------initilization of N particles--------------------
	initials=zeros(9,n);
	for i=1:n
		phi = rand*2*pi;   % yaw
        chi = rand*2*pi;     % pitch
        psi = rand*2*pi;     % roll
        
		pose=zeros(3,3);
		pose(1,1)=cos(phi)*cos(chi);
		pose(1,2)=cos(phi)*sin(chi)*sin(psi)-sin(phi)*cos(psi);
		pose(1,3)=cos(phi)*sin(chi)*cos(psi)+sin(phi)*sin(psi);
		pose(2,1)=sin(phi)*cos(chi);
		pose(2,2)=sin(phi)*sin(chi)*sin(psi)+cos(phi)*cos(psi);
		pose(2,3)=sin(phi)*sin(chi)*cos(psi)-cos(phi)*sin(psi);
		pose(3,1)=-sin(chi);
		pose(3,2)=cos(chi)*sin(psi);
		pose(3,3)=cos(chi)*cos(psi);
 
		initials(:,i)=reshape(pose,9,1);
	end;

% 	initials=rand(9,n)*2-1;
%    for i=1:n
%		i_initial=initials(:,i);
%		i_matrix=reshape(i_initial,3,3);
%	    [U,sigma,V]=svd(i_matrix);
%		i_rotation=U*V';
%		initials(:,i)=reshape(i_rotation,9,1);
%	end;
	current_states=initials;
       	
	%-------------------------------------------------------
	k=0;	
	t=1;
    weights=(1/n)*ones(1,n);   

	while t<T
		k=k+1
		t=k^power;                       % temperature at each iteration
        %global_pro=1/(1+k^0.1);          % probability of using global proposal
		global_pro=0.5;
		delta_t_k=k^1.1-(k-1)^1.1;       % temperature gap can be different at different iteration
     		
        max_objective=-inf; 
	    for i=1:n	
		    i_current_state=current_states(:,i);	
			i_current_rotation=reshape(i_current_state,3,3);        
			i_current_objective=computing_objective(s_source,target,i_current_rotation);	
			
			% coputes the importance weights
			i_ratio=(i_current_objective+delta)^delta_t_k;
			weights(i)=weights(i)*i_ratio;

			if i_current_objective>max_objective
				max_objective=i_current_objective;
				optimal_rotation=i_current_rotation;
			end;
		end; 

		% compute the relative effective sample size
		weights=weights/sum(weights);
		%square_of_sum=(sum(weights))^2;
		sum_of_square=sum(weights.^2);
	    ESS=1/(n*sum_of_square)
		%save('current_states.mat','current_states');
	 
		% computer the diversity of particles
	    std_deviation=std(current_states'); 	
		mean_deviation=mean(std_deviation)
	    	

		% if ESS is smaller than a threshold or c 
		if (ESS<ESS_threshold) || (mean_deviation<0.5)
			% disp('in resampling and MCMC transition')
			% systematic resampling 
			if ESS<ESS_threshold
				current_states=systematic_resample(current_states,weights,n); 
			end;
				
			weights=(1/n)*ones(1,n);
			% MCMC transition 
			for i=1:n
                r=rand;
				if r<global_pro
					i_next_state=2*rand(9,1)-1;
				else
					i_current_state=current_states(:,i);	
			    	i_current_rotation=reshape(i_current_state,3,3);        
			    	i_current_objective=computing_objective(s_source,target,i_current_rotation);		
					i_next_state=mvnrnd(i_current_state, (nu^2)*eye(9))+mvnrnd(zeros(9,1), gamma^2*eye(9));
				end;

				i_next_matrix=reshape(i_next_state,3,3);
				[U,sigma,V]=svd(i_next_matrix);
				i_next_rotation=U*V';
				i_next_objective=computing_objective(s_source, target, i_next_rotation);

	        	prob=((i_next_objective+delta)/(i_current_objective+delta))^t;
				
				r=rand;
				if r<prob 
	%				disp('accept');
					i_next_state=reshape(i_next_rotation,9,1);
					current_states(:,i)=i_next_state;
				end;	
			end;
		end; 

		max_objective
		trace=[trace,max_objective];
	end; 
elseif strcmp(method,'ASMC')
    %--------------initilization of N particles--------------------
    %initials=zeros(9,n);
	%for i=1:n
	%	phi = rand*2*pi;   % yaw
    %    chi = rand*pi;     % pitch
    %    psi = rand*pi;     % roll
    %    
	%	pose=zeros(3,3);
	%	pose(1,1)=cos(phi)*cos(chi);
	%	pose(1,2)=cos(phi)*sin(chi)*sin(psi)-sin(phi)*cos(psi);
	%	pose(1,3)=cos(phi)*sin(chi)*cos(psi)+sin(phi)*sin(psi);
	%	pose(2,1)=sin(phi)*cos(chi);
	%	pose(2,2)=sin(phi)*sin(chi)*sin(psi)+cos(phi)*cos(psi);
	%	pose(2,3)=sin(phi)*sin(chi)*cos(psi)-cos(phi)*sin(psi);
	%	pose(3,1)=-sin(chi);
	%	pose(3,2)=cos(chi)*sin(psi);
	%	pose(3,3)=cos(chi)*cos(psi);
 
	%	initials(:,i)=reshape(pose,9,1);
	%end;

	initials=rand(9,n)*2-1;
    for i=1:n
		i_initial=initials(:,i);
		i_matrix=reshape(i_initial,3,3);
	    [U,sigma,V]=svd(i_matrix);
		i_rotation=U*V';
		initials(:,i)=reshape(i_rotation,9,1);
	end;
	current_states=initials;
	current_objectives=zeros(1,n);
       	
	%-------------------------------------------------------
	k=0;	
	t=1;
	num_history=floor(size_subsample/n);
	sub_chain=cell(1,num_history);   % 10 previous chain states are recoreded for adapting covaraince 

	while t<T
		k=k+1
		
		pro_adaptive=1/k^(0.5);    % the probability of adapting covariance 

		weights=zeros(1,n);   
        max_objective=-inf; 
	    for i=1:n	
		    i_current_state=current_states(:,i);	
			i_current_rotation=reshape(i_current_state,3,3);        
			i_current_objective=computing_objective(s_source,target,i_current_rotation);	
			current_objectives(i)=i_current_objective;
	
			
			if i_current_objective>max_objective
				max_objective=i_current_objective;
				optimal_rotation=i_current_rotation;
			end;
		end; 
		idx=mod(k,num_history);
		if idx==0
			idx=num_history;
		end;
		sub_chain{1,idx}=current_states;

		% compute the delta_t_k for the next iteration
	    delta_t_k=find_optimal_gap(current_objectives, delta, t, T, ESS_threshold)
		t=t+delta_t_k;
        
		% computes the importance weights
		weights(i)=(i_current_objective+delta)^delta_t_k;	

		% systematic resampling 
        current_states=systematic_resample(current_states,weights,n); 

	    % MCMC transition 
		% adapt the covriance 
        if k>=10
			r=rand; 
            r
            pro_adaptive
			if (r<pro_adaptive)
				history_data=zeros(9,num_history*n);
				for tt=1:num_history
					history_data(:,(tt-1)*n+1:(tt-1)*n+n)=sub_chain{1,tt};
				end;
                empirical_covariance=cov(history_data');
			    mean_diag=mean(diag(empirical_covariance));
				empirical_covariance=empirical_covariance/mean_diag;
			end
			%	history_data=zeros(9,k*n);
			%	for tt=1:k
			%		histtory_data(:,(tt-1)*n+1:(tt-1)*n+n)=sub_chain{1,tt};
			%	end;
		else
			empirical_covariance=eye(9);
		end;

	    for i=1:n
           i_current_state=current_states(:,i);	
	       i_current_rotation=reshape(i_current_state,3,3);        
	       i_current_objective=computing_objective(s_source,target,i_current_rotation);	

		   i_next_state=mvnrnd(i_current_state, (gamma^2)*eye(9)+(nu^2)*empirical_covariance);
	   	   i_next_matrix=reshape(i_next_state,3,3);
	   	   [U,sigma,V]=svd(i_next_matrix);
	   	   i_next_rotation=U*V';
	   	   i_next_objective=computing_objective(s_source, target, i_next_rotation);

	   	   prob=((i_next_objective+delta)/(i_current_objective+delta))^t;
	   	   
	   	   r=rand;
	   	   if r<prob 	
		%	   disp('accepted');
			   i_next_state=reshape(i_next_rotation,9,1);
	   	   	   current_states(:,i)=i_next_state;
	   	   end;	
	    end;
	    

	    max_objective
	    trace=[trace,max_objective];
	end; 
elseif (strcmp(method,'KASMC-g'))
   
%elseif (method=='KASMC')
end;






%********************************************************************************************
function [objective]=computing_objective(s_source,target,rotation)
num_point=size(s_source,1);

s_transformed=rotation*s_source';
s_transformed=s_transformed';

% ----search the nearest neighbour in target for each point in s_transformed
idx=knnsearch(target,s_transformed);     % idx contains the index of nearest neighbour in target for the corresponding row in s_transformed
idx=idx(1:num_point); 
nearest_neighbours=target(idx,:);

error=0;
for i=1:num_point
	error=error+norm(s_transformed(i,:)-nearest_neighbours(i,:));
end;
%error=norm(nearest_neighbours-s_transformed,'fro');

error=error/num_point;

objective=1/(1+error);




%*******************************************************************************************
function [r_current_states]=systematic_resample(current_states,weights,num)
[dim,num_data]=size(current_states);
% construct output data
r_current_states=zeros(dim,num_data);

v=num_data*weights;
j=1;
c=v(1);
u=rand;

% at first resample the num_data samples
for k=1:num_data
	while c<u
		if j+1 >num_data
			break
		end;
		j=j+1;
		c=c+v(j);
	end;

	r_current_states(:,k)=current_states(:,j);

	u=u+1;
end;

% extract a subset (with num size) of resampled data
if num<num_data
	rand_idx=randperm(num_data);
	r_current_states=r_current_states(:,rand_idx(1:num));
end;
	

function [r_current_states]=multinomial_resample(current_states,weights,num)
[dim,num_data]=size(current_states);
% construct output data
r_current_states=zeros(dim,num);

weights=weights/(sum(weights));
cum_weights=cumsum(weights);

for i=1:num
	r=rand;
	idx=find(cum_weights>=r,1);
	r_current_states(:,i)=current_states(:,idx);
end;


function [delta_t_k]=find_optimal_gap(objectives, delta, J_k, J, ESS_threshold)
num_data=size(objectives,2);
l=0;
u=50;
a=(l+u)/2;
while 1		
	weights=(objectives+delta).^a;
	% compute the relative effective sample size  (ESS)
	square_of_sum=(sum(weights))^2;
	sum_of_square=sum(weights.^2);
	ESS=square_of_sum/(num_data*sum_of_square);

	if ESS<ESS_threshold
		u=a;
        a=(l+u)/2;
	else
		l=a;
        a=(l+u)/2;
	end;

	if (abs(u-l)<0.005)
		delta_t_k=min(a,J-J_k);
    	break;
	end;
end;


