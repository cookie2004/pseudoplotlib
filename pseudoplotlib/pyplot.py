# import matplotlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects
from matplotlib import animation
from matplotlib.gridspec import GridSpec 
from matplotlib import collections as mcoll

from pseudoplotlib.utils import *
from pseudoplotlib.colormaps import *

def kabsch(a, b, weights=None, return_v=False):
  '''Get alignment matrix for two sets of coordinates, with optional weights.'''
  
  if weights is not None:
    # Apply weights to the coordinates
    weights = np.sqrt(weights)[:, np.newaxis]
    a_weighted = a * weights
    b_weighted = b * weights
  else:
    a_weighted = a
    b_weighted = b
  
  # Compute the weighted covariance matrix
  ab = a_weighted.swapaxes(-1, -2) @ b_weighted
  
  # Perform Singular Value Decomposition
  u, s, vh = np.linalg.svd(ab, full_matrices=False)
  
  # Handle reflection case
  flip = np.linalg.det(u @ vh) < 0
  u_ = np.where(flip, -u[...,-1].T, u[...,-1].T).T
  u[...,-1] = u_
  
  return u if return_v else (u @ vh)

def nankabsch(a,b,**kwargs):
    """
    Get alignment matrix for two sets of coordinates, excluding rows with NaN or infinite values.
    
    Parameters:
    a (ndarray): First set of coordinates.
    b (ndarray): Second set of coordinates.
    **kwargs: Additional keyword arguments passed to the kabsch function.
    
    Returns:
    ndarray: Alignment matrix.
    """
    # Check for finite values in both sets of coordinates
    ok = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    # Filter out rows with NaN or infinite values
    a,b = a[ok],b[ok]
    # Call the kabsch function with the filtered coordinates and additional arguments
    return kabsch(a,b,**kwargs)

def plot_pseudo_3D(xyz, c=None, ax=None, chainbreak=5, Ls=None,
                   cmap="gist_rainbow", line_w=2.0,
                   cmin=None, cmax=None, zmin=None, zmax=None,
                   shadow=0.95,remove_axes=False):
  """
  Plots a pseudo-3D representation of a protein structure in 2D, with optional color mapping.

  Parameters:
  xyz (array-like): Coordinates of the protein structure.
  c (array-like, optional): Colors for each segment.
  ax (matplotlib.axes.Axes, optional): Matplotlib Axes object to plot on.
  chainbreak (float, optional): Distance to consider a break in the chain.
  Ls (list of int, optional): Lengths of chains.
  cmap (str or Colormap, optional): Colormap to use for coloring segments.
  line_w (float, optional): Line width.
  cmin, cmax (float, optional): Min and max values for color scaling.
  zmin, zmax (float, optional): Min and max values for z-dimension scaling.
  shadow (float, optional): Shadow intensity.

  Returns:
  matplotlib.collections.LineCollection: Line collection added to the plot.
  """
  # make segments and colors for each segment
  xyz = np.asarray(xyz)
  if Ls is None:
    seg = np.concatenate([xyz[:,None],np.roll(xyz,1,0)[:,None]],axis=1)
    c_seg = np.arange(len(seg))[::-1] if c is None else (c + np.roll(c,1,0))/2
  else:
    Ln = 0
    seg = []
    c_seg = []
    for L in Ls:
      sub_xyz = xyz[Ln:Ln+L]
      seg.append(np.concatenate([sub_xyz[:,None],np.roll(sub_xyz,1,0)[:,None]],axis=1))
      if c is not None:
        sub_c = c[Ln:Ln+L]
        c_seg.append((sub_c + np.roll(sub_c,1,0))/2)
      Ln += L
    seg = np.concatenate(seg,0)
    c_seg = np.arange(len(seg))[::-1] if c is None else np.concatenate(c_seg,0)
  
  # set colors
  c_seg = rescale(c_seg,cmin,cmax)  
  if isinstance(cmap, str):
    if cmap == "gist_rainbow": 
      c_seg *= 0.75
    colors = matplotlib.cm.get_cmap(cmap)(c_seg)
  else:
    colors = cmap(c_seg)
  
  # remove segments that aren't connected
  seg_len = np.sqrt(np.square(seg[:,0] - seg[:,1]).sum(-1))
  if chainbreak is not None:
    idx = seg_len < chainbreak
    seg = seg[idx]
    seg_len = seg_len[idx]
    colors = colors[idx]

  seg_mid = seg.mean(1)
  seg_xy = seg[...,:2]
  seg_z = seg[...,2].mean(-1)
  order = seg_z.argsort()

  # add shade/tint based on z-dimension
  z = rescale(seg_z,zmin,zmax)[:,None]

  # add shadow (make lines darker if they are behind other lines)
  seg_len_cutoff = (seg_len[:,None] + seg_len[None,:]) / 2
  seg_mid_z = seg_mid[:,2]
  seg_mid_dist = np.sqrt(np.square(seg_mid[:,None] - seg_mid[None,:]).sum(-1))
  shadow_mask = sigmoid(seg_len_cutoff * 2.0 - seg_mid_dist) * (seg_mid_z[:,None] < seg_mid_z[None,:])
  np.fill_diagonal(shadow_mask,0.0)
  shadow_mask = shadow ** shadow_mask.sum(-1,keepdims=True)

  seg_mid_xz = seg_mid[:,:2]
  seg_mid_xydist = np.sqrt(np.square(seg_mid_xz[:,None] - seg_mid_xz[None,:]).sum(-1))
  tint_mask = sigmoid(seg_len_cutoff/2 - seg_mid_xydist) * (seg_mid_z[:,None] < seg_mid_z[None,:])
  np.fill_diagonal(tint_mask,0.0)
  tint_mask = 1 - tint_mask.max(-1,keepdims=True)

  colors[:,:3] = colors[:,:3] + (1 - colors[:,:3]) * (0.50 * z + 0.50 * tint_mask) / 3
  colors[:,:3] = colors[:,:3] * (0.20 + 0.25 * z + 0.55 * shadow_mask)

  set_lim = False
  if ax is None:
    fig, ax = plt.subplots()
    fig.set_figwidth(5)
    fig.set_figheight(5)
    set_lim = True
  else:
    fig = ax.get_figure()
    if ax.get_xlim() == (0,1):
      set_lim = True
      
  if set_lim:
    xy_min = xyz[:,:2].min() - line_w
    xy_max = xyz[:,:2].max() + line_w
    ax.set_xlim(xy_min,xy_max)
    ax.set_ylim(xy_min,xy_max)

  ax.set_aspect('equal')
    
  # determine linewidths
  width = fig.bbox_inches.width * ax.get_position().width
  linewidths = line_w * 72 * width / np.diff(ax.get_xlim())

  lines = mcoll.LineCollection(seg_xy[order], colors=colors[order], linewidths=linewidths,
                               path_effects=[matplotlib.patheffects.Stroke(capstyle="round")])

  if remove_axes:
      ax.axis('off')
  return ax.add_collection(lines)

def plot_ticks(ax, Ls, Ln=None, add_yticks=False):
  if Ln is None: Ln = sum(Ls)
  L_prev = 0
  for L_i in Ls[:-1]:
    L = L_prev + L_i
    L_prev += L_i
    ax.plot([0,Ln],[L,L],color="black")
    ax.plot([L,L],[0,Ln],color="black")
  
  if add_yticks:
    ticks = np.cumsum([0]+Ls)
    ticks = (ticks[1:] + ticks[:-1])/2
    ax.yticks(ticks,alphabet_list[:len(ticks)])

def make_animation(xyz,
                   seq=None,
                   sitewise=None,
                   pairwise=None,
                   
                   sitewise_label="plddt",
                   sitewise_min=0.5,
                   sitewise_max=0.9,
                   sitewise_color=None,
                   
                   pairwise_label="pae",
                   pairwise_min=0.0,
                   pairwise_max=30.0,
                   
                   losses=None,
                   pos_ref=None,
                   line_w=2.0,
                   dpi=100, interval=60, color_msa="Taylor",
                   length=None, align_xyz=True, **kwargs):
  """
  Creates an animation of protein folding trajectories.

  Parameters:
  xyz (list of np.ndarray): List of 3D coordinates for each frame.
  seq (list of np.ndarray): Sequence data for each frame.
  sitewise (list of np.ndarray): Sitewise data (e.g., confidence scores) for each frame.
  pairwise (list of np.ndarray): Pairwise data (e.g., distance matrices) for each frame.
  sitewise_label (str): Label for the sitewise data.
  sitewise_min (float): Minimum value for sitewise data color scaling.
  sitewise_max (float): Maximum value for sitewise data color scaling.
  sitewise_color (str): Color scheme for sitewise data.
  pairwise_label (str): Label for the pairwise data.
  pairwise_min (float): Minimum value for pairwise data color scaling.
  pairwise_max (float): Maximum value for pairwise data color scaling.
  losses (list of float): Loss values for each frame.
  pos_ref (np.ndarray): Reference position for alignment.
  line_w (float): Line width for plotting.
  dpi (int): Dots per inch for the figure.
  interval (int): Interval between frames in milliseconds.
  color_msa (str): Color scheme for multiple sequence alignment.
  length (int or list): Length of the protein or list of lengths for each chain.
  align_xyz (bool): Whether to align the coordinates to a reference.
  kwargs: Additional keyword arguments.

  Returns:
  matplotlib.animation.ArtistAnimation: The generated animation.
  """                     
  if pos_ref is None:
    pos_ref = xyz[-1]

  if length is None:
    L = len(pos_ref)
    Ls = None
  elif isinstance(length, list):
    L = length[0]
    Ls = length
  else:
    L = length
    Ls = None

  # align to reference
  if align_xyz:
      
    pos_ref_trim = pos_ref[:L]
    pos_ref_trim_mu = np.nanmean(pos_ref_trim,0)
    pos_ref_trim = pos_ref_trim - pos_ref_trim_mu

    # align to reference position
    new_pos = []
    for x in xyz:
      x_mu = np.nanmean(x[:L],0)
      aln = nankabsch(x[:L]-x_mu, pos_ref_trim)
      new_pos.append((x-x_mu) @ aln)

    pos = np.array(new_pos)

    # rotate for best view
    pos_mean = np.concatenate(pos,0)
    m = np.nanmean(pos_mean,0)
    rot_mtx = nankabsch(pos_mean - m, pos_mean - m, return_v=True)
    pos = (pos - m) @ rot_mtx
    pos_ref_full = ((pos_ref - pos_ref_trim_mu) - m) @ rot_mtx
  
  else:
    # rotate for best view
    pos_mean = np.concatenate(xyz,0)
    m = np.nanmean(pos_mean,0)
    aln = nankabsch(pos_mean - m, pos_mean - m, return_v=True)
    pos = [(x - m) @ aln for x in xyz]
    pos_ref_full = (pos_ref - m) @ aln

  # initialize figure
  fig = plt.figure()
  gs = GridSpec(4,3, figure=fig)
  if pairwise is None:
    if seq is None:
      ax1 = fig.add_subplot(gs[:,:])
    else:
      ax1, ax2 = fig.add_subplot(gs[:3,:]), fig.add_subplot(gs[3:,:])
  else:
    if seq is None:
      ax1, ax3 = fig.add_subplot(gs[:,:2]), fig.add_subplot(gs[:,2:])
    else:
      ax1, ax2, ax3 = fig.add_subplot(gs[:3,:2]), fig.add_subplot(gs[3:,:]), fig.add_subplot(gs[:3,2:])

  fig.subplots_adjust(top=0.95,bottom=0.1,right=0.95,left=0.05,hspace=0,wspace=0)
  fig.set_figwidth(8); fig.set_figheight(6); fig.set_dpi(dpi)
  if seq is not None:
    if seq[0].shape[0] > seq[0].shape[1]:
      ax2.set_ylabel("positions")
      ax2.set_xlabel("sequences")
      ax2.set_xticks([])
    else:
      ax2.set_xlabel("positions")
      ax2.set_ylabel("sequences" if seq[0].shape[0] > 1 else "amino acids")
      ax2.set_yticks([])

  if sitewise is None:
    ax1.set_title("N→C")
  else:
    ax1.set_title(sitewise_label)
  
  if pairwise is not None:
    ax3.set_title(pairwise_label)
    ax3.set_xticks([])
    ax3.set_yticks([])

  # set bounderies
  main_pos = pos_ref_full[np.isfinite(pos_ref_full).all(1)]
  pred_pos = [np.isfinite(x).all(1) for x in pos]
  x_min,y_min,z_min = np.minimum(np.mean([x.min(0) for x in pred_pos],0),main_pos.min(0)) - 5
  x_max,y_max,z_max = np.maximum(np.mean([x.max(0) for x in pred_pos],0),main_pos.max(0)) + 5

  x_pad = ((y_max - y_min) * 2 - (x_max - x_min)) / 2
  y_pad = ((x_max - x_min) / 2 - (y_max - y_min)) / 2
  if x_pad > 0:
    x_min -= x_pad
    x_max += x_pad
  else:
    y_min -= y_pad
    y_max += y_pad

  ax1.set_xlim(x_min, x_max)
  ax1.set_ylim(y_min, y_max)
  ax1.set_xticks([])
  ax1.set_yticks([])

  # get animation frames
  ims = []
  for k in range(len(pos)):
    ims.append([])
    flags = dict(ax=ax1, line_w=line_w, zmin=z_min, zmax=z_max)
    if sitewise is None:
      if sitewise_color == "chain":
        c = np.concatenate([[n]*L for n,L in enumerate(length)])
        ims[-1].append(plot_pseudo_3D(pos[k], c=c,  Ls=Ls, cmap=pymol_cmap, cmin=0, cmax=39, **flags))
      else:
        L = pos[k].shape[0]
        ims[-1].append(plot_pseudo_3D(pos[k], c=np.arange(L)[::-1],  Ls=Ls, cmin=0, cmax=L, **flags))  
    else:
      ims[-1].append(plot_pseudo_3D(pos[k], c=sitewise[k], Ls=Ls, cmin=sitewise_min, cmax=sitewise_max, **flags))

    if seq is not None:
      if seq[k].shape[0] == 1:
        ims[-1].append(ax2.imshow(seq[k][0].T, animated=True, cmap="bwr_r",vmin=-1, vmax=1))
      else:
        cmap = matplotlib.colors.ListedColormap(jalview_color_list[color_msa])
        vmax = len(jalview_color_list[color_msa]) - 1
        msa_oh = seq[k][:,:,:20]
        msa = msa_oh.argmax(-1).astype(float)
        msa[msa_oh.sum(-1) == 0] = np.nan
        if msa.shape[0] > msa.shape[1]:
          msa = msa.T
        ims[-1].append(ax2.imshow(msa, animated=True, cmap=cmap, vmin=0, vmax=vmax, interpolation="none"))
    
    if pairwise is not None:
      L = pairwise[k].shape[0]
      ims[-1].append(ax3.imshow(pairwise[k], animated=True, cmap="bwr",
                                vmin=pairwise_min, vmax=pairwise_max, extent=(0, L, L, 0)))

  # add lines
  if length is not None:
    Ls = length if isinstance(length, list) else [length,None]
    if pairwise is not None:
      plot_ticks(ax3, Ls, pairwise[0].shape[0])

  # make animation!
  ani = animation.ArtistAnimation(fig, ims, blit=True, interval=interval)
  plt.close()
  return ani
