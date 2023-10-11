function [] = damk()

Das = [0.3, 2e5];
colors = [0 0.4470 0.7410; 
    0.4660, 0.6740, 0.1880];

R = 8.314; % J/mol/K
T = 50:50:2000; % C
rho0 = 2780; % kg/m3
t0 = 50e6 * 3.154e7; % seconds

T0 = 2000; % C

Sj = 0.5 .* 0.2 .* 0.1 .* 0.2; % 0.001
dg = 0.001; % grain size, m

n = 1; % mol
rj = exp(-T0./T);

for i = 1:numel(Das)
    Da = Das(i);
    col = colors(i,:);
    Gamma0 = Da .* rho0 / t0; % kg/m3/s

    % min Aj
    Aj = 80; % J
    Rj = Gamma0 .* rj .* Sj .* Aj ./ n ./ R ./ T; % kg/m3/s
    Rj_surface = Rj ./ dg; % kg/m2/s

    Rj_gcy = Rj_surface .*1e3 ./ 100 ./100 .* 3.154e7; % g/cm2/yr

    % max Aj
    
    Aj = 5600; % J 
    Rj = Gamma0  .* rj .* Sj .* Aj ./ n ./ R ./ T;
    Rj_surface = Rj ./ dg;
    Rj_gcy2 = Rj_surface .*1e3 ./ 100 ./100 .* 3.154e7; % g/cm2/yr
    pshape = polyshape([T flip(T)], [log10(Rj_gcy) flip(log10(Rj_gcy2))]);
    p=plot(pshape, 'EdgeColor','none','FaceAlpha',0.3,'FaceColor',col);
    hold on

    % med Aj
    
    Aj = 675; % J 
    Rj = Gamma0 .* rj .* Sj .* Aj ./ n ./ R ./ T;
    Rj_surface = Rj ./ dg;
    Rj_gcy3 = Rj_surface .*1e3 ./ 100 ./100 .* 3.154e7; % g/cm2/yr
    plot(T, log10(Rj_gcy3),'--','Color',col,'LineWidth',1.5);
end

ylim([-11,1])
xlim([350,950])
yticks(-11:2:0)

saveas(gcf,'da.pdf')

end