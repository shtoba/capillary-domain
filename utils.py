#!/usr/bin/env python
# coding: utf-8

# In[77]:


get_ipython().run_line_magic('matplotlib', 'inline')

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pandas as pd
import random
import collections
import pprint
import sys
from io import StringIO
from statistics import stdev, mean


# In[1]:


# Second ver. 200124

def points_to_lines(points):
    all_lines = []

    field_borders = [{'id': 'bottom', 'a': 0, 'b': 1, 'c': -field_min_y},
                     {'id': 'top', 'a': 0, 'b': 1, 'c': -field_max_y},
                     {'id': 'left', 'a': 1, 'b': 0, 'c': -field_min_x},
                     {'id': 'right', 'a': 1, 'b': 0, 'c': -field_max_x}]


    for i, p1 in enumerate(points):
        lines_tmp = []

        other_points = points.copy()
        other_points.remove(p1)

        for p2 in other_points:
            dict_tmp = {'id': [i, points.index(p2)]}

            a = (p2[0] - p1[0])
            b = (p2[1] - p1[1])
            c = - (p2[0]**2 - p1[0]**2 + p2[1]**2 - p1[1]**2) / 2

            dict_tmp['a'] = a
            dict_tmp['b'] = b
            dict_tmp['c'] = c

            lines_tmp.append(dict_tmp)

        lines_tmp.extend(field_borders)
        all_lines.append(lines_tmp)
    
    return all_lines


# In[79]:


def get_intersection(line1, line2):
    if line1['a'] * line2['b'] == line2['a'] * line1['b']:
        return
    else:
        x = (
            ((line1['b'] * line2['c']) - (line2['b'] * line1['c']))
            / ((line1['a'] * line2['b']) - (line2['a'] * line1['b']))
        )

        y = (
            ((line2['a'] * line1['c']) - (line1['a'] * line2['c']))
            / ((line1['a'] * line2['b']) - (line2['a'] * line1['b']))
        )

        dict_tmp = {'x': x, 'y': y}
    
        return dict_tmp


# In[80]:


def get_y(line, x):
    if line['b'] == 0:
        return
    else:
        return -(line['a'] * x + line['c']) / line['b']


# In[2]:


def lines_to_intersections(all_lines):
    all_intersections = []
    
    # for each point
    for k, lines in enumerate(all_lines):
        intersections_tmp = []
        lines_done = []
        
        # for each line
        for line1 in lines:    
            other_lines = lines.copy()
            
            lines_done.append(line1)
            for l_done in lines_done:
                other_lines.remove(l_done)
            
            for line2 in other_lines:
                intxn = get_intersection(line1, line2)
                if intxn:
                    dict_tmp = {'id': [line1['id'], line2['id']]}
                    dict_tmp.update(intxn)
                    intersections_tmp.append(dict_tmp)
        
        all_intersections.append(intersections_tmp)
            
    return all_intersections


# In[11]:


def get_nearest_intersections(points, all_lines, all_intersections, show=False, show_count=True):
    if show_count: print('*** get_nearest_intersections() ***')
    nearest_intersections = []
    
    for point_id, intersections in enumerate(all_intersections):
        if show:
            print('*** Point id:', point_id)
        else:
            if show_count: print('\rPoint id: {} / {}'.format(point_id+1, len(points)), end='')
            
        p = points[point_id]
        nearest_intersections_tmp = []
        
        for intxn in intersections:
            if show: print('    Intxn:', intxn)
            line_to_intxn = {'a': intxn['y'] - p[1],
                             'b': p[0] - intxn['x'],
                             'c': (intxn['x'] - p[0]) * p[1] - (intxn['y'] - p[1]) * p[0]}
            
            lines = [l for l in all_lines[point_id] if not l['id'] in intxn['id']]
            
            add_to_list = True
            
            for line in lines:
                if show: print('        Line:', line)
                intxn_new = get_intersection(line, line_to_intxn)
                if show: print('          intxn_new:', intxn_new, end='')
                
                if intxn_new:
                    min_x = min(p[0], intxn['x'])
                    max_x = max(p[0], intxn['x'])
                    min_y = min(p[1], intxn['y'])
                    max_y = max(p[1], intxn['y'])

                    if min_x <= intxn_new['x'] <= max_x and min_y <= intxn_new['y'] <= max_y:
                        if intxn_new['y'] == intxn['y'] and intxn_new['x'] == intxn['x']:
                            if show: print(' On the line')
                        else:
                            if show: print(' NG')
                            add_to_list = False
                    else:
                        if show: print('')
                else:
                    if show: print('')
            
            if add_to_list:
                existing_intxn_list = [[a['x'], a['y']] for a in nearest_intersections_tmp]
                if [intxn['x'], intxn['y']] in existing_intxn_list:
                    if show: print('        -> NG (Already exists)', end='\n\n')
                else:
                    if show: print('        -> OK', end='\n\n')
                    nearest_intersections_tmp.append(intxn)
            else:
                if show: print('        -> NG (Not nearest)', end='\n\n')
    
        nearest_intersections.append(nearest_intersections_tmp)
        if show: pprint.pprint(nearest_intersections_tmp)
        if show: print('')

    return nearest_intersections


# In[83]:


def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]


# In[4]:


def get_line_segments(nearest_intersections, show=False):
    if show: print('')
    if show: print('*** get_line_segments() ***')

    segments = []
    
    nearest_intersections = organize_intersections(nearest_intersections)
    
    for point_id, intxns in enumerate(nearest_intersections):
        if show: print('')
        if show: print('*** Point id:', point_id)
        segments_tmp = []
        
        intxns.append(intxns[0])
        
        for i in range(len(intxns) - 1):
            x1 = intxns[i]['x']
            x2 = intxns[i+1]['x']
            y1 = intxns[i]['y']
            y2 = intxns[i+1]['y']
            
            if show: print('    segment:({}, {}) - ({}, {})'.format(x1, y1, x2, y2))
            
            segments_tmp.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
            
        segments.append(segments_tmp)
    
    return segments


# In[5]:


def organize_intersections(nearest_intersections):
    intxns_new = []
        
    for point_id, intxns in enumerate(nearest_intersections):
        intxns.sort(key=lambda i: i['x'])
        first = intxns[0]
        last = intxns[-1]
        
        line = {'a': first['y'] - last['y'],
                'b': last['x'] - first['x'],
                'c': (first['x'] - last['x']) * last['y'] - (first['y'] - last['y']) * last['x']}
        
        unders = []
        overs = []
        
        saved_coords = []
        
        for intx in intxns:
            if not [intx['x'], intx['y']] in saved_coords:
                saved_coords.append([intx['x'], intx['y']])
                if intx['y'] > get_y(line, intx['x']):
                    overs.append(intx)
                else:
                    unders.append(intx)
        
        unders = unders[::-1]
        overs.extend(unders)
        
        intxns_new.append(overs)
                
    return intxns_new


# In[6]:


def plot(points, all_lines, nearest_intersections, segments, simple=False, fill=False):
    plt.xlim(field_min_x, field_max_x)
    plt.ylim(field_min_y, field_max_y)
    
    coord_x = [x[0] for x in points]
    coord_y = [x[1] for x in points]
    
    colors = cm.gist_rainbow(np.linspace(0, 1, len(points)))
    
    delta_ratio = 0.02
    delta_x = np.array([(x[0]-(field_max_x-field_min_x)/2)*delta_ratio for x in points])
    delta_y = np.array([(x[1]-(field_max_y-field_min_y)/2)*delta_ratio for x in points])
    
    if not simple:
        for point_id, lines in enumerate(all_lines):
            for line in lines:
                if line['b'] == 0:
                    plt.plot([-line['c'] / line['a'], -line['c'] / line['a']],
                             [field_min_y, field_max_y],
                             c='tab:gray',
                             linewidth=0.1)
                else:
                    plt.plot([field_min_x, field_max_x],
                             [get_y(line, field_min_x), get_y(line, field_max_x)],
                             c='tab:gray',
                             linewidth=0.1)

        
    if not simple:
        for point_id, intxn in enumerate(nearest_intersections):
            #rnd_value = np.array([random_value[point_id] for i in range(len(intxn))])
            rnd_value_x = np.array([delta_x[point_id] for i in range(len(intxn))])
            rnd_value_y = np.array([delta_y[point_id] for i in range(len(intxn))])
            coord_x_intxn = np.array([x['x'] for x in intxn])
            coord_y_intxn = np.array([x['y'] for x in intxn])
            plt.scatter(coord_x_intxn + rnd_value_x, coord_y_intxn + rnd_value_y,
                        c=[colors[point_id]],
                        #alpha=1,
                        edgecolors=[colors[point_id]],
                        clip_on=False,
                        s=50)
    
    for point_id, segs in enumerate(segments):
        if point_id == 7 or True:
            for seg in segs:
                if simple:
                    plt.plot([seg['x1'],
                              seg['x2']],
                             [seg['y1'],
                              seg['y2']],
                             c=colors[point_id], clip_on=False, alpha=0.6, linewidth=2)
                    
                else:
                    plt.plot([seg['x1'] + delta_x[point_id],
                              seg['x2'] + delta_x[point_id]],
                             [seg['y1'] + delta_y[point_id],
                              seg['y2'] + delta_y[point_id]],
                             c=colors[point_id], clip_on=False, alpha=0.6, linewidth=2)
                    
    if fill:
        nearest_intersections = organize_intersections(nearest_intersections)
        for point_id, intxns in enumerate(nearest_intersections):
            coords = [[i['x'], i['y']] for i in intxns]
            tmp = plt.Polygon(coords, fc=colors[point_id], alpha=0.15)
            plt.gca().add_patch(tmp)
                
    for i in range(len(points)):
        plt.scatter(coord_x[i], coord_y[i], c=[colors[i]], marker='o', s=100, label=i)     
        
    plt.legend(bbox_to_anchor=[1.2,1], ncol=len(points)//13+1)
        
    plt.show()
    plt.close()


# In[7]:


def calculate_areas(nearest_intersections):
    areas_list = []
    
    nearest_intersections = organize_intersections(nearest_intersections)
    
    for point_id, intxns in enumerate(nearest_intersections):
        x = [i['x'] for i in intxns]
        y = [i['y'] for i in intxns]
        
        area = 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
        areas_list.append(area)
        
    field = (field_max_x - field_min_x) * (field_max_y - field_min_y) / 1000 / 1000
    cap_dens = len(nearest_intersections) / field
    
    areas = np.array(areas_list)
    radius = np.sqrt(areas / np.pi)
    radius_sd = stdev(radius)
    
    radius_log = np.log10(radius)
    radius_log_sd = stdev(radius_log)
    
    plt.hist(radius)
    plt.xlabel('Radius')
    plt.ylabel('Frequency')
    plt.gca().get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    plt.show()
    plt.close()
    
    print('Density (No./mm2)    :', cap_dens)
    print('Mean radius (um)     :', mean(radius))
    print('SD radius (um)       :', radius_sd)
    print('Mean log radius (um) :', mean(radius_log))
    print('SD log radius (x1000):', radius_log_sd * 1000)

    results = np.stack([areas, radius, radius_log], 1)
    df = pd.DataFrame(results, columns=['area', 'radius', 'radius_log'])
    
    return df


# In[8]:


def get_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


# In[ ]:




