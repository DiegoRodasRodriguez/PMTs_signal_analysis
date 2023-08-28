function Y = RamseyRecoModel(Qinf, Qo, K, E, titleLabel)

Y = Qo + (Qinf-Qo)./(1+K./E);

if nargin>4
    plot(Y, E, '-');
    title(titleLabel);   
end

end