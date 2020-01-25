#!/usr/bin/env python
# coding: utf-8

# In[60]:


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


# In[2]:


# Second ver. 200124

def points_to_lines(points):
    lines = {}

    field_borders = [{'id': 'bottom', 'a': 0, 'b': 1, 'c': -field_min_y},
                     {'id': 'top', 'a': 0, 'b': 1, 'c': -field_max_y},
                     {'id': 'left', 'a': 1, 'b': 0, 'c': -field_min_x},
                     {'id': 'right', 'a': 1, 'b': 0, 'c': -field_max_x}]


    for i, p1 in enumerate(points):
        lines[i] = []

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

            lines[i].append(dict_tmp)

        lines[i].extend(field_borders)
    
    return lines


# In[3]:


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


# In[4]:


def get_y(line, x):
    if line['b'] == 0:
        return
    else:
        return -(line['a'] * x + line['c']) / line['b']


# In[5]:


def lines_to_intersections(all_lines):
    intersections = {}
    
    # for each point
    for k, lines in all_lines.items():
        intersections[k] = []
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
                    intersections[k].append(dict_tmp)
            
    return intersections


# In[6]:


def get_nearest_intersections(points, all_lines, all_intersections, show=False):
    nearest_intersections = {}
    
    for point_id, intersections in all_intersections.items():
        if show: print('*** Point id:', point_id)
        p = points[point_id]
        nearest_intersections[point_id] = []
        duplicate = False
        
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
                existing_intxn_list = [[a['x'], a['y']] for a in nearest_intersections[point_id]]
                if [intxn['x'], intxn['y']] in existing_intxn_list:
                    duplicate = True
                    if show: print('        -> OK (Already exists)', end='\n\n')
                else:
                    if show: print('        -> OK', end='\n\n')

                nearest_intersections[point_id].append(intxn)
            else:
                if show: print('        -> NG', end='\n\n')
    
        if show: pprint.pprint(nearest_intersections[point_id])
        if show: print('')

        # Algorithm for removing duplicated intersections (may not be perfect)
        if duplicate:
            line_ids = [x['id'][0] for x in nearest_intersections[point_id]]
            line_ids.extend([x['id'][1] for x in nearest_intersections[point_id]])
            if show: print('line_ids', line_ids)
            count = collections.Counter(tuple(item) for item in line_ids)
            
            duplicate_ids_org = []
            for k, v in count.most_common():
                if v > 2:
                    for i in range(v - 2):
                        duplicate_ids_org.append(k)
            
            if show: print(count)
            if show: print('duplicate_ids_org', duplicate_ids_org)
            
            finish = False
            
            while not finish:
                if show: print('')
                if show: print('*** New loop ***')
                
                duplicate_ids = random.sample(duplicate_ids_org, len(duplicate_ids_org))
                
                if show: print('duplicate_ids', duplicate_ids)
                
                n_intxns_copy = nearest_intersections[point_id].copy()
                
                for id1 in duplicate_ids:
                    if show: print('id1', id1)
                    other_duplicate_ids = duplicate_ids.copy()
                    other_duplicate_ids = [x for x in other_duplicate_ids if not x == id1]

                    for id2 in other_duplicate_ids:
                        if show: print('id2', id2)
                        if len(id1) > 2:
                            id1_org = ''.join(id1)
                        else:
                            id1_org = list(id1)

                        if len(id2) > 2:
                            id2_org = ''.join(id2)
                        else:
                            id2_org = list(id2)

                        intxn_to_remove = [d for d in nearest_intersections[point_id]                                           if d['id'] == [id1_org, id2_org] or d['id'] == [id2_org, id1_org]]
                        
                        if show: print(intxn_to_remove)

                        if intxn_to_remove:
                            try:
                                duplicate_ids.remove(id1)
                                duplicate_ids.remove(id2)

                                if show: print('Removing:', intxn_to_remove)

                                n_intxns_copy.remove(intxn_to_remove[0])

                                if show: print('After removal')
                                if show: pprint.pprint(n_intxns_copy)

                                if not duplicate_ids:
                                    finish = True
                                    if show: print('')
                                    if show: print('Done')
                                    
                                    nearest_intersections[point_id] = n_intxns_copy
                                    
                                    if show: pprint.pprint(nearest_intersections[point_id])
                                    if show: print('')

                            except:
                                break
                    else:
                        continue
                    break
        
    return nearest_intersections


# In[7]:


def get_unique_list(seq):
    seen = []
    return [x for x in seq if x not in seen and not seen.append(x)]


# In[8]:


def get_line_segments(nearest_intersections):
    segments = {}
    
    for point_id, intxns in nearest_intersections.items():
        segments[point_id] = []
        
        line_ids = [x['id'][0] for x in intxns]
        line_ids.extend([x['id'][1] for x in intxns])
        line_ids = get_unique_list(line_ids)
        
        for line_id in line_ids:
            #print(line_id)
            xs = [i['x'] for i in intxns if line_id in i['id']]
            ys = [i['y'] for i in intxns if line_id in i['id']]
            
            #if line_id == [1, 3]:
            #    print('xs', xs)
            #    print('ys', ys)
                
            id_1 = xs.index(min(xs))
            id_2 = 0 if id_1 == 1 else 1 
            
            #print(xs)
            #print(id_2)
            
            x1 = xs[id_1]
            y1 = ys[id_1]
            x2 = xs[id_2]
            y2 = ys[id_2]
            
            segments[point_id].append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
    
    return segments


# In[74]:


def organize_intersections(nearest_intersections):
    intxns_new = {}
        
    for point_id, intxns in nearest_intersections.items():
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
        
        intxns_new[point_id] = overs
                
    return intxns_new


# In[67]:


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
        for point_id, lines in all_lines.items():
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
        for point_id, intxn in nearest_intersections.items():
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
    
    for point_id, segs in segments.items():
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
        for point_id, intxns in nearest_intersections.items():
            coords = [[i['x'], i['y']] for i in intxns]
            tmp = plt.Polygon(coords, fc=colors[point_id], alpha=0.15)
            plt.gca().add_patch(tmp)
                
    for i in range(len(points)):
        plt.scatter(coord_x[i], coord_y[i], c=[colors[i]], marker='o', s=100, label=i)     
        
    plt.legend(bbox_to_anchor=[1.2,1], ncol=len(points)//10+1)
        
    plt.show()
    plt.close()


# In[75]:


def calculate_areas(nearest_intersections):
    areas_list = []
    
    nearest_intersections = organize_intersections(nearest_intersections)
    
    for point_id, intxns in nearest_intersections.items():
        x = [i['x'] for i in intxns]
        y = [i['y'] for i in intxns]
        
        area = 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
        areas_list.append(area)
        
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
    
    print('Mean radius     :', mean(radius))
    print('SD radius       :', radius_sd)
    print('Mean log radius :', mean(radius_log))
    print('SD log radius   :', radius_log_sd)

    results = np.stack([areas, radius, radius_log], 1)
    df = pd.DataFrame(results, columns=['area', 'radius', 'radius_log'])
    
    return df


# In[15]:


def get_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


# In[ ]:




