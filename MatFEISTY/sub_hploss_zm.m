%%% Offline coupling
function [zf] = sub_hploss_zm(enc_1,enc_2,enc_3,enc_4,enc_5,bio_1,bio_2,bio_3,bio_4,bio_5,dZ)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    con_1 = enc_1 .* bio_1;
    con_2 = enc_2 .* bio_2;
    con_3 = enc_3 .* bio_3;
    con_4 = enc_4 .* bio_4;
    con_5 = enc_5 .* bio_5;
    
    % Fraction of zooplankton mortality loss consumed
    zf = (con_1 + con_2 + con_3 + con_4 + con_5) ./ (dZ+eps);

end
