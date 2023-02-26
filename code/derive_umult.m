function N = derive_umult(N,P)
    if P.umult==0 %estimates based on range of ice thicknesses
        [N.umult,N.sigUmult,~]=f_probability2(N.THX(N.MASK)); %optimizes based on range of ice thicknesses
%         umult=0.9; %range is [0.8,1], more realistically [0.81,0.99]
    elseif umult==2 %estimates based on ice thickness of each individual pixel
        [N.umult,N.sigUmult]=f_probability1(N.THX,N.MASK); %optimizes based on each ice thickness
    elseif umult<0.8|umult>1 
        error('Unphyiscal multiplier for column-averaged velocity')
    else
        N.umult=P.umult;
        N.sigUmult=0.1;
    end
