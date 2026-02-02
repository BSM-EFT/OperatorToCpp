import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm

label_dict = {
    'mut3': r'$m_{\tilde{t}_R}$',
    'm1': r'$m_1$'
    }

## Functions to convert a given WC name into LaTeX

def texify(wc_name):
    parts = wc_name.split("_")
    if len(parts) == 2:
        [name, idx] = parts
        if name[-1] >= '1' and name[-1] <= '9':
            idx = "{" + '(' + name[-1] + ')' + idx + "}" 
            name = "{" + name[1:-1] + "}"
        else:
            name = "{" + name[1:] + "}"
            idx = "{" + idx + "}"
        return [name, idx]
    else:
        name = parts[0][1:]
        if name == "HBox":
            name = "{" + "H\square" + "}"
        else:
            name = "{" + name + "}"
        return [name]

def texify_abs(wc_name):   
    parts = texify(wc_name)
    if len(parts) == 2:
        [name, idx] = parts
        tex = r'$|C_%s^%s\,\, ({\rm TeV}^{-2})|$'%(name,idx)
    else:
        tex = r'$|C_%s\,\, ({\rm TeV}^{-2})|$'%(parts[0])
    
    return tex

def texify_pos(wc_name):   
    parts = texify(wc_name)
    if len(parts) == 2:
        [name, idx] = parts
        tex = r'$C_%s^%s\,\, ({\rm TeV}^{-2})$'%(name,idx)
    else:
        tex = r'$C_%s\,\, ({\rm TeV}^{-2})$'%(parts[0])
    
    return tex

def texify_neg(wc_name):   
    parts = texify(wc_name)
    if len(parts) == 2:
        [name, idx] = parts
        tex = r'$-C_%s^%s\,\, ({\rm TeV}^{-2})$'%(name,idx)
    else:
        tex = r'$-C_%s\,\, ({\rm TeV}^{-2})$'%(parts[0])
    
    return tex

def texify_nounits(wc_name):   
    parts = texify(wc_name)
    if len(parts) == 2:
        [name, idx] = parts
        tex = r'$|C_%s^%s|$'%(name,idx)
    else:
        tex = r'$|C_%s|$'%(parts[0])
    
    return tex

def tex_label(data, wc_name):
    is_neg = [num < 0.0 for num in np.array(data[wc_name].array) ]
    is_pos = [num >= 0.0 for num in np.array(data[wc_name].array) ]
    if all(is_neg):
        return texify_neg(wc_name)
    elif all(is_pos):
        return texify_pos(wc_name)
    else:
        return texify_abs(wc_name)
    
# function for making 2d plots

def plot_fn_2d(data, param_names, wc_name, filename):
    assert len(param_names) == 2

    param_vals = []
    for p_name in param_names:
        param_vals.append(np.unique(data[p_name].array))
    param_vals = np.array(param_vals)

    x_axis = param_vals[0]
    y_axis = param_vals[1]

    x_label = label_dict[param_names[0]] + r'$\,\, {\rm (TeV)}$'
    y_label = label_dict[param_names[1]] + r'$\,\, {\rm (TeV)}$'

    cbar_label = tex_label(data, wc_name)

    wc_vals = abs(np.array(data[wc_name].array))
    vals = wc_vals.reshape(len(x_axis),len(y_axis))

    fig, ax = plt.subplots(figsize=(8, 6))
    img = ax.imshow(vals, cmap='afmhot', norm=LogNorm(),
                    origin='lower', extent=[x_axis.min(), x_axis.max(), y_axis.min(), y_axis.max()])

    ax.plot([x_axis.min(), x_axis.max()], [y_axis.min(), y_axis.max()], color='gray', linestyle='--', linewidth=1.1)
    x_fill = np.linspace(x_axis.min(), x_axis.max(), 100)
    ax.fill_between(x_fill, x_fill, y_axis.max(), color='#AAAAAA', alpha=0.4)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(cbar_label, fontsize=27)
    ax.set_xlabel(x_label, fontsize=27)
    ax.set_ylabel(y_label, fontsize=27)
    ax.set_title(r'$g_{1} = 0.37, \quad g_{3} = 1.1$', fontsize=27, pad=16)
    ax.set_xticks([0.5,1.0,1.5,2.0,2.5])
    ax.set_yticks([0.5,1.0,1.5,2.0,2.5])
    plt.tight_layout()
    plt.savefig(filename)

def collect_vals(data, bmk_pts, wc_names):
    values = []
    for bp in bmk_pts:
        cond = True
        for (k,v) in bp.items():
            cond = cond & (data[k] == v)

        d1 = data[cond]
        
        for i, row in d1.iterrows():
            values.append(abs(row[wc_names].values))
    
    return values

def create_legend(bp):
    cmds = ""
    exprs = []
    for k in bp.keys():
        cmds += "%s $= {\\rm %s\\, TeV}$, "
        exprs.append(label_dict[k])
        exprs.append(bp[k])

    cmds = cmds[:-2]
    exprs = tuple(exprs)
    legend = cmds%exprs
    
    return legend

# function for making bar plots

def plot_fn_bar(data, bmk_pts, wc_names, filename,**kwargs):
    x_pos = np.arange(len(wc_names))
    fig, ax = plt.subplots(figsize=(15, 8))

    bar_width = kwargs.get('bar_width')
    y_min = kwargs.get('y_min')
    y_max = kwargs.get('y_max')

    vals = collect_vals(data, bmk_pts, wc_names)

    ax.bar(x_pos - bar_width, vals[2], width=bar_width, label=create_legend(bmk_pts[2]), color='#FFFFAA', edgecolor='black', alpha=0.85)
    ax.bar(x_pos, vals[1], width=bar_width, label=create_legend(bmk_pts[1]), color='tab:red', edgecolor='black', alpha=0.75)
    ax.bar(x_pos + bar_width, vals[0], width=bar_width, label=create_legend(bmk_pts[0]), color='#030764', edgecolor='black', alpha=0.85)

    ax.set_yscale('log')
    ax.set_ylim(y_min,y_max)
    ax.set_xticks(x_pos)

    ticklabels = []
    for wc in wc_names:
        ticklabels.append(texify_nounits(wc))

    ax.set_xticklabels(ticklabels)  
    ax.set_ylabel(r'${\rm TeV}^{-2}$', fontsize=24)
    ax.legend(fontsize=20,ncol=1,framealpha=1.0,loc=4)

    plt.tight_layout()
    plt.savefig(filename)
