function [A,A_plus,A_minus,R_plus,R_minus,C,p]=BA(m0,m,N,p,mu,sigma)
    myOnes=ones(m0,1);
    i=[1:m0];
    j=[(m0+1)*ones(1,m0)];
    mk=N;
    nk=N;
    s=ones(m0,1);
    %The total nb of connections is 2*m*NbIterations
    Net=sparse(i,j,s,mk,nk,N*m);
     %Vector contains the degree of nodes
    D=zeros(N,1);
    %Connects the m0+1'th vertex
    D(1:m0)=myOnes;
    D(m0+1)=m;

    %buffers,temp identifies the choosed vertices


    %Incremental growth of the network,now the network has m0+i-1 vertices
    for i=m0+2:N
        Prob=cumsum(D/sum(D));
        randm=rand(1,m);
        temp=zeros(1,i-1);
        for j=1:m%choose m vertices 
            for k=1:i-1%from the exsist i-1 vertices
                if(randm(j)<=Prob(k)&temp(k)==0)
                    Net(k,i)=1;
                    D(k)=D(k)+1;
                    D(i)=D(i)+1;
                    temp(k)=1;
                    break;
                end
            end
        end
    end

    %full
    Net=full(Net);
%     Net=Net+Net';
%     D=D';
    nonzero_indices = find(Net ~= 0);
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
    s=1;
    A=R_plus-R_minus-s*diag(ones(1,N));
    C=(sum(A(:)~=0)-N)/(N*(N-1));
end