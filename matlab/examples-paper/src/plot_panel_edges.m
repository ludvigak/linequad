function plot_panel_edges(g, h, varargin)

if isempty(varargin)
    style = '-k';
else
    style = varargin{1};
end

for i=1:g.numpanels+1
    plot(complex(g.z_edges(i) + g.n_edges(i)*[h -h]), style)
    hold on
end
