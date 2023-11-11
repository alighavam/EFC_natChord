function c = map2color(vec, cmap)

% Normalize your values to the colormap range
c = interp1(linspace(min(vec), max(vec), size(cmap, 1)), cmap, vec);

