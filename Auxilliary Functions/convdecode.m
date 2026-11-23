function [r]=convdecode(x,k)

    if k==2 %1/2
    
        G=[56721 61713];
        %G=[133 171];
        %G=[4335 5723];
        %G=[561 753];
        %G=[447254 627324]
        %G=[716502 514576]
        L=15;
    
    elseif k==3 %R=2/3
    
        G=[63 15 46;
        32 65 61];
    
        L=[6 6];
    
    elseif k==4 %r=3/4
    
        G=[40 14 34 60;
        04 64 20 70;
        34 00 60 64];
    
        L=[6 6 6];
    end

    trellis=poly2trellis(L, G);

    %tb = ceil(length(x)/2);  
    tb=10*max(L);

    r=vitdec(x, trellis, tb, 'trunc', 'hard');

end

