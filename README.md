# Overview
This code is truing (straightening) narrow curved tube-like structure, determined by spot positions (but in principle those can be pixels). Idea is to split this structure on to several segments, and each segment wit to the polynomial curve of second order. Segments should be overlapped with each other for stitching afterwards.

The usage is quite complex, so if someone wants to try ut I recommend to contact me first for initial run together.

# Usage 
## Step 1 
Firstly we need isolate only spots which belongs to CNS. for this I use annoations and annotation list of ROIs which are part of CNS. This is usually done with notebook: **get_whole_cns2.ipynb**

## Step 2
Run curve truing for only a segment of CNS. for this chose segment that will be well fitted with one parabola (so only one curvatuer point). Also please move this curvature point in such way that parabola will look either upside or downside (not laterally). The most important are parameteres in configuration file (see **conf_truing.yaml**):
 - *xlim*, *ylim*: limits for the segment of the CNS. if the tube segment is complex and cant be simply fitted in rectangular, you can take larger rectangular and then crop out small rectangulars in *subcrop_x*, *subcrop_y*
 - *pos_a*, *pos_b*: approximate initial positions of outline spots to start outline determination
 - *vector_a*, *vector_b*: initial vectors drawn from *pos_a*, *pos_b* to determine outline of the CNS
 - out_folder: output folder
 - rotate_90: whether or not rotate on 90 degrees whole segment
 - shorten_median_prc: on how much (in %) shorten the median determination (from the outline determination), as closer to the end there's more probability of wrong median determination. Check computed median from output images and check if it well represents center of the tube
 - method: can be "default" or "skeleton". This is the method of median determination (this is the most complicated part of the code). I usually use default, skeleton give often many branches
 - rotate_angle: custom rotational angle for whoel curve to make it looks like **U** or inverted **U**
 - rotate_c_point: only if rotation angle is not 0, then I use center of selected region (*xlim*, *ylim*)
As output it produces 2 images (curve_plot and fit plot), where I usually check quality of: median building, and fitting this median to a polyfit curve

## Step 3
When you straighten first two segments with as in **Step 2**, you will need to stitch them horizontally using **stitching_sections.py** code. Here are also important parameters specified in configuration file (**conf_stitching.yaml**):
 - *path_section_right_csv*, *path_section_left_csv*: paths to output csvs from truing step, make sure you place them correctly (right_csv should be on the right hand side after stitching and left at the left)
 - *flip_x_left*, *flip_y_left*: as "left" is usually a new section, you can flip it horizintally or vertically before stitching if needed
 - *center_search_shift*: both sections have some overlap. this is very approximately their overlap value in position units (usually pixels).
 - *range_search_shift*: program will go +- this value from *center_search_shift* while searching for the best shift value. THe best shift value is simply where the distance between the same spots in both segments is minimal. Please pay attention that final shift should be in the range *center_search_shift* +- *range_search_shift*, not at one of the ends of this range.
 - *path_output*: where csv file with merged spots positions will be saved

# General procedure
For one CNS tube you have to start with **Step 1**, and then sequentially perform **Step 2** with **Step 3**: which means truing new segment, and stitch it to all previously truined and merged segments. You have to check performance of the code at each step and correct parameters in configuration files. I think in average I re-run truing curve code for the same segmen 5-10 times with different parameteres before I become fully satisfied with the results

# TODO
1) Add environment file
2) Structural code change - can we do it all at once, without breaking it to segments?
