function plot_solution(q,f,x,t)
    figure
    subplot(211)
    contourf(t,x,real(q),'DisplayName','q','linecolor','none');
    colorbar();
    xlabel('t');
    ylabel('x');
    legend();
    subplot(212)
    plot(x,sum(abs(f).^2,2));
    xlabel('t');
    ylabel('Force Magnitude');
end