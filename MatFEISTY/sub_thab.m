%%% Habitat Temp
function temp = sub_thab(Tp,Tb,tdif)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));

end