def plot_velocity(self, **kwargs):
    raise NotImplementedError("Coming soon!")
    start_timestep = kwargs.pop('start_timestep', 0)
    stop_timestep  = kwargs.pop('stop_timestep', len(self.time))
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    ax = axes.quiver(self.x, self.y, self.u[0,:], self.v[:,0], vmin=vmin, vmax=vmax, **kwargs)
    axes.axis('scaled')
    def update(i):
       ax.set_array(self.Dataset['zs'][i,:-1,:-1].ravel())
       return ax
    
    anim = FuncAnimation(fig, update, frames=np.arange(start_timestep+1, stop_timestep), interval=interval)
    if colorbar==True:
        plt.colorbar(ax)
    if show is True:
        plt.show()
    return anim