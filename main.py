#!/usr/bin/env python
# coding: utf-8

# # Initialization

# In[43]:


get_ipython().run_line_magic('run', 'utils.ipynb')


# In[19]:


field_min_x = 0
field_max_x = 15
field_min_y = 0
field_max_y = 10

points0 = [[1,5],
           [2,2],
           [2,8],
           [4,8],
           [5,1],
           [6,6],
           [10,4],
           [10,6],
           [10,8],
           [12,9],
           [12,1],
           [14,6]]

points1 = [[2,8],
          [10,8],
          [6,6]]

points2 = [[2,8],
          [10,8],
          [6,2],
          [6,4],
          [6,6]]

points3 = [[4,3],[4,6],
           [8,3],[8,6],
           [12,3],[12,6]]

points4 = [[1,1],
          [1,4],
          [4,1],
          [4,4],
          [2,2],
          [3,2]]


# In[20]:


import random
rnd_x = sorted([random.random() * 15 for i in range(50)])
rnd_y = [random.random() * 10 for i in range(50)]
points5 = [[x, y] for x, y in zip(rnd_x, rnd_y)]


# # Calculation

# In[21]:


points = points3

all_lines = points_to_lines(points)
all_intersections = lines_to_intersections(all_lines)
nearest_intersections = get_nearest_intersections(points,
                                                  all_lines,
                                                  all_intersections,
                                                  show=False)
segments = get_line_segments(nearest_intersections)


# # Results

# In[44]:


plot(points, all_lines, nearest_intersections, segments)
plot(points, all_lines, nearest_intersections, segments, simple=True, fill=True)


# In[45]:


organize_intersections(nearest_intersections)


# In[46]:


#sys.stdout = stdout_backup # print log for get_nearest_intersections()
areas = calculate_areas(nearest_intersections)
areas


# In[ ]:




