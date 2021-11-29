d50 = load ('data_60uni_1meV.mat');
d50_uni = load ('data_60uni_5meV.mat');
d50_mix = load ('data_60uni_mix.mat');

d50_pss = d50.P_ss_all;
d50_uni_pss = d50_uni.P_ss_all;
d50_mix_pss = d50_mix.P_ss_all;

d50_Len = d50.Len;
d50_uni_Len = d50_uni.Len;
d50_mix_Len= d50_mix.Len;

figure(1)
plot([2:60], d50_uni_pss(2:60), 'x-b')
hold
plot([2:60], d50_pss(2:60), 'x-r')
plot([2:60], d50_mix_pss(2:60), 'x-g')
legend({'uniform','mix'},'Location','northeast')
ylabel('Flux','interpreter','latex')
xlabel('N')
newcolors = {'#0aa303','#3189ff'};
colororder(newcolors)
% set(gca, 'YScale', 'log')
set(gca, 'FontName', 'Times New Roman','FontSize', 18)