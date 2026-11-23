function [y]=digqammod(x,N)
    y=[];
    if (N==4| N==16| N==64| N==256) & mod(length(x),log2(N))==0
        M=log2(N);
        for i=1:M:length(x)-M+1
            y=[y,bi2de(x(i:i+M-1),'left-msb')];
        end
        y=qammod(y,N,'Gray');
    end
end


