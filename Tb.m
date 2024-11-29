function [PX,PY,T_values]=Tb(xi,C,JZ,FC,s)

    T = @(z, xi, C, JZ, FC,s) FC*sum((xi.^2)./((xi+z).*(xi+conj(z)))) - (s+JZ)^2;
%      T = @(z, xi, C, JZ, FC,s) FC*sum((xi.^2)./((xi+z+JZ).*(xi+JZ+conj(z)))) - (s+JZ)^2;
    deta_x=max(xi)-min(xi);
    rl=-max(xi)-deta_x;
    rr=-min(xi)+deta_x;
    
    
%     x = linspace(-3, 1, 2000);
    x = linspace(rl, rr, 1000);
    y = linspace(-0.1, 0.1, 101);
    
    [X, Y] = meshgrid(x, y);
    Z = X + 1i * Y;
    T_values = arrayfun(@(z) T(z, xi, C, JZ, FC,s), Z);
    
    result = false(size(T_values));
    
    for row = 1:size(T_values, 1)
        % 找到第一个正值的索引
        cf = find(T_values(row, :) > 0);
        if numel(cf)>0
            result(row,cf(1))=true;
            result(row,cf(end))=true;
        end
    end
    for col = 1:size(T_values, 2)
        % 找到第一个正值的索引
        cf = find(T_values(:, col) > 0);
        if numel(cf)>0
            result(cf(1),col)=true;
            result(cf(end),col)=true;
        end
    end
    PX=X(result)*(s+JZ);
    PY=Y(result)*(s+JZ);

end
