function [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s)

temp=rand(N,N)<C;
temp=triu(temp,1);
nonzero_indices = find(temp ~= 0);
num_elements = round(p * numel(nonzero_indices));
selected_indices = randsample(nonzero_indices, num_elements);
remaining_indices = setdiff(nonzero_indices, selected_indices);
A_plus=zeros(N,N);
A_minus=zeros(N,N);
A_plus(selected_indices)=1;
A_plus=A_plus+A_plus';
A_minus(remaining_indices)=1;
A_minus=A_minus+A_minus';
R=sigma.*randn(N,N)+mu;
R_plus=A_plus.*R;
R_minus=A_minus.*R;
A=R_plus-R_minus-s*eye(N); 

end