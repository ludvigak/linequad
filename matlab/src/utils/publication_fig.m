function publication_fig
set(gcf,'PaperPositionMode','auto')
P = get(gcf,'Position');
set(gcf,'Position',[P(1) P(2) 275 220])
set(gca,'FontName','Times','FontSize',8)
box off
