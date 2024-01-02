function [E_an,i_ter] = kepler_e_9(M_an,ecc_1)
    E_an0 = M_an;
    E_an1 = M_an + ecc_1*sin(E_an0);
    i_ter = 1;
    E_an = M_an + ecc_1*sin(E_an1);
    while ( abs(E_an-E_an1) >= 10^-9)
       E_an1 = E_an;
       E_an = M_an + ecc_1*sin(E_an1);
       i_ter = i_ter +1;
    end
