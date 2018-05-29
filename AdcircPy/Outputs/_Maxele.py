def make_plot(self, extent=None, axes=None, title=None, step=0, total_colors=256, cbar_label=None, **kwargs):
    axes, idx = _plotters._init_fig(self, axes, extent, title)
    if isinstance(self.values, list):
        values = self.values[step]
    else:
        values = self.values
    vmin = kwargs.pop("vmin", np.min(values))
    vmax = kwargs.pop("vmax", np.max(values))
    cmap = kwargs.pop("cmap", "jet")
    levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
    if np.ma.is_masked(values):
        trimask = np.any(values.mask[self.elements], axis=1)
        Tri = matplotlib.tri.Triangulation(self.x, self.y, self.elements, trimask)
        axes.tricontourf(Tri, values, levels=levels, cmap=cmap, extend='both')
    else:
        axes.tricontourf(self.x, self.y, self.elements, values, levels=levels, cmap=cmap, extend='both')
    cbar = _plotters._init_colorbar(axes, cmap, vmin, vmax)
    if cbar_label is not None:
        cbar.set_label(cbar_label)
    cbar.set_ticks([vmin,
                    vmin+(1./4.)*(vmax-vmin),
                    vmin+(1./2.)*(vmax-vmin),
                    vmin+(3./4.)*(vmax-vmin),
                    vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 
                            np.around(vmin+(1./4.)*(vmax-vmin), 2),
                            np.around(vmin+(1./2.)*(vmax-vmin), 2),
                            np.around(vmin+(3./4.)*(vmax-vmin), 2),
                            np.around(vmax, 2)])
    return axes