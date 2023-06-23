from scipy import ndimage
import numpy as np

def zoomcontour(ax,xgrid,ygrid,zgrid,**kwargs):
    xz = ndimage.zoom(xgrid, 3.0)
    yz = ndimage.zoom(ygrid, 3.0)
    zz = ndimage.zoom(zgrid, 3.0)
    
    return ax.contour(xz, yz, zz, **kwargs)
