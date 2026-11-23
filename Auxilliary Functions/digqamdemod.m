function [x]=digqamdemod(y,N)
    if (N==4| N==16| N==64| N==256)
        x=qamdemod(y,N,'Gray');
        x=de2bi(x,log2(N),'left-msb');
        x=reshape(x',1,[]);
    end
end


