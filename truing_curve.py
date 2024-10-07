import scanpy as sc
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fire
import yaml
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from skimage.morphology import skeletonize


def reorder_polygon_points(pos_points, start_n = 0):
    start_n = np.where(pos_points[:,0]==np.max(pos_points[:,0]))
    start_n = start_n[0][0]
    print(start_n)
    current_point = pos_points[start_n]
    pos_points = np.delete(pos_points, start_n, 0)
    new_set_points = []; new_set_points.append(tuple(current_point))
    prev_point = current_point
    while pos_points.shape[0]>2:
        current_point, n_cur = find_closest_point(prev_point, pos_points, False)
        new_set_points.append(tuple(current_point))
        pos_points = np.delete(pos_points, n_cur, 0)
        prev_point = current_point
    return new_set_points




def find_middle_point(point_a, point_b):
    middle_point = np.array([(point_a[0]+point_b[0])/2, (point_a[1]+point_b[1])/2])
    return middle_point

def find_next_point(prev_point, vector, set_of_points):
    proj_point = prev_point+vector
    next_point, n = find_closest_point(proj_point, set_of_points, False)
    next_vector = next_point-prev_point
    return next_point, next_vector, n

def iterate_through_points(point_a, vector_a, point_b, vector_b, all_pos, first_median_point, perc_to_cut):
    median_points = []; median_points.append(first_median_point)
    N_steps_approx = int(all_pos.shape[0]/2 - all_pos.shape[0]*perc_to_cut/100)
    for i in range(N_steps_approx):
        point_a,vector_a, n_a = find_next_point(point_a,vector_a, all_pos)
        all_pos = np.delete(all_pos, n_a, 0)
        point_b,vector_b, n_b = find_next_point(point_b,vector_b, all_pos)
        
        #if i<10:
        print(i)
        print(point_a)
        print(point_b)
        print(vector_a)
        print(vector_b)
        
        all_pos = np.delete(all_pos, n_b, 0)
        median_points.append(find_middle_point(point_a, point_b))
        #print('New point a is ' + str(point_a))
        #print('New point b is ' + str(point_b))
    return np.array(median_points), point_a,vector_a,point_b,vector_b, all_pos

def find_first_2_points(start_point, set_of_points):
    point_a, n_a = find_closest_point(start_point, set_of_points)
    vector_a = point_a - start_point
    set_of_points2 = np.delete(set_of_points, n_a, 0)
    point_b, n_b = find_closest_point(start_point, set_of_points2)
    vector_b = point_b - start_point
    return point_a, point_b, vector_a, vector_b, n_a, n_b


def find_parabola_xc(x0, A, B, C):
    #C1 = C+A*x0**2-B*x0
    #B1 = B-2*A
    #x_c = (-B1 - np.sqrt(B1**2-4*A*C1))/(2*A)
    x_c = -(B-2*x0*A)/(2*A)
    return x_c

def get_arc_length(x_list, y_list, n1, n2):
    Length = 0
    diff_x = np.diff(x_list); diff_y = np.diff(y_list)
    if n1>n2:
        n3 = n2; n2=n1; n1 = n3
    for i in range(n1,n2):
        Length += np.sqrt(1+(diff_y[i]/diff_x[i])**2)* np.abs(diff_x[i])
    return Length


def find_closest_point(point, set_of_points, include0):
    d = []
    for i in range(set_of_points.shape[0]): d.append(math.dist(point,set_of_points[i]))
    d = np.array(d); 
    if not include0: d[d==0] = 999999999999 #just to make sure we dont pickup the same point
    n_min = np.argmin(d)
    pos_min = set_of_points[n_min]
    return pos_min, n_min

def straighten_curve(pts_x, pts_y, fit_x, fit_y, n_c):
    str_pts_set = np.zeros((pts_x.shape[0],2))
    set_fit_pts = np.zeros((fit_x.shape[0],2)); set_fit_pts[:,0] = fit_x; set_fit_pts[:,1] = fit_y
    x_c = fit_x[n_c]; i=0
    for x,y in zip(pts_x,pts_y):
        fit_closest_pt, n = find_closest_point(np.array([x,y]), set_fit_pts, True)
        
        dist = np.sqrt((x-fit_closest_pt[0])**2 + (y-fit_closest_pt[1])**2)
        #here i check whether point is "under" the parabola or above
        point = Point(np.array([x,y])); polygon = Polygon(set_fit_pts)
        if not polygon.contains(point):
            dist*=(-1)
            
        arc_length = get_arc_length(fit_x, fit_y, n_c, n)
        if x>x_c:
            str_pts_set[i] = np.array([arc_length+x_c, dist])
        else:
            str_pts_set[i] = np.array([x_c-arc_length, dist])
        i+=1
        #print(i, end='\r')
    return str_pts_set

def get_roi_outline(path_adata):
    adata=sc.read_h5ad(path_adata)
    ROI = adata.obs[adata.obs['cns']==1]
    outline = adata.obs[adata.obs['outline']==1]
    return ROI, outline

def ReadConfFile(FilePath):
    with open(FilePath, 'r') as file:
        data = yaml.safe_load(file)
    path_adata = data['path_adata']
    name = data['name']
    x_lim = data['x_lim']
    y_lim = data['y_lim']
    pos_a = np.array(data['pos_a'])
    pos_b = np.array(data['pos_b'])
    v_a = np.array(data['vector_a'])
    v_b = np.array(data['vector_b'])
    out_folder = data['out_folder']
    rotate_90 = data['rotate_90']
    U_shape = data['U_shape']
    save_images = data['save_images']
    shorten_median_prc = data['shorten_median_prc']
    subcrop_x = data['subcrop_x']
    subcrop_y = data['subcrop_y']
    downscale = data['downscale']
    method = data['method']
    rotate_angle = data['rotate_angle']
    rotate_c_point = data['rotate_c_point']
    return path_adata, name, x_lim, y_lim, pos_a, pos_b, out_folder, rotate_90, U_shape, save_images, shorten_median_prc, subcrop_x, subcrop_y, downscale, method, rotate_angle, rotate_c_point, v_a, v_b

def crop_roi_outline(ROI, outline, x_lim, y_lim):
    ROI_part = ROI[ROI['spatial_1']>x_lim[0]]; ROI_part = ROI_part[ROI_part['spatial_1']<x_lim[1]]
    ROI_part = ROI_part[ROI_part['spatial_2']>y_lim[0]]; ROI_part = ROI_part[ROI_part['spatial_2']<y_lim[1]]
    outline_part = outline[outline['spatial_1']>x_lim[0]]; outline_part = outline_part[outline_part['spatial_1']<x_lim[1]]
    outline_part = outline_part[outline_part['spatial_2']>y_lim[0]]; outline_part = outline_part[outline_part['spatial_2']<y_lim[1]]
    return ROI_part, outline_part

def subcrop_roi_outline(ROI, outline, crop_x, crop_y):

    i=0
    while i<ROI.shape[0]:
        if ROI['spatial_1'][i]>crop_x[0] and ROI['spatial_1'][i]<crop_x[1]:
            if ROI['spatial_2'][i]>crop_y[0] and ROI['spatial_2'][i]<crop_y[1]:
                #print(ROI['spatial_1'][i])
                #print(ROI['spatial_2'][i])
                ROI.drop(ROI.index[i], inplace=True)
            else: 
                i+=1
        else:
            i+=1
    i=0
    while i<outline.shape[0]:
        if outline['spatial_1'][i]>crop_x[0] and outline['spatial_1'][i]<crop_x[1]:
            if outline['spatial_2'][i]>crop_y[0] and outline['spatial_2'][i]<crop_y[1]:
                outline.drop(outline.index[i], inplace=True)
            else:
                i+=1
        else:
            i+=1        
    
    #ROI_crop = ROI[ROI['spatial_1']<crop_x[0]]; ROI_crop = ROI_crop[ROI_crop['spatial_1']>crop_x[1]]
    #ROI_crop = ROI_crop[ROI_crop['spatial_2']<crop_y[0]]; ROI_crop = ROI_crop[ROI_crop['spatial_2']>crop_y[1]]
    #outline_crop = outline[outline['spatial_1']<crop_x[0]]; outline_crop = outline_crop[outline_crop['spatial_1']>crop_x[1]]
    #outline_crop = outline_crop[outline_crop['spatial_2']<crop_y[0]]; outline_crop = outline_crop[outline_crop['spatial_2']>crop_y[1]]
    return ROI, outline



def make_binary_image_from_outline(pos_outline, downscale = 100):
    Xmax = np.max(pos_outline[:,0]); Ymax = np.max(pos_outline[:,1]); Max = np.max([Xmax, Ymax])
    img = np.zeros((int(Max/downscale), int(Max/downscale)), dtype = 'uint8')
    #n_start = 
    polygon_points = reorder_polygon_points(pos_outline/downscale)
    polygon = Polygon(polygon_points)
    
    for i in range(img.shape[0]):
        print(i, end = '\r')
        for j in range(img.shape[1]):
            point = Point(i, j)
            if polygon.contains(point):
                img[i,j] = 1
    return img

def define_median_modes(all_pos_outline, pos_a, pos_b, v_a, v_b, shorten_median_prc, out_folder, name, downscale, mode = 'skeleton'):
    #print('v_a = ' + str(v_a))
    
    if mode=='skeleton':
        #downscale = 200
        img = make_binary_image_from_outline(all_pos_outline, downscale)
        skeleton = skeletonize(img, method = 'lee')
        
        
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                         sharex=True, sharey=True)

        ax = axes.ravel()

        ax[0].imshow(img, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('original', fontsize=20)

        ax[1].imshow(skeleton, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('skeleton', fontsize=20)
        path_fig = out_folder + '/' + name + '_binary_img.png'
        plt.savefig(path_fig, format="jpg", dpi = 150)
        
        y_m, x_m = np.nonzero(skeleton)
        median_pts = np.zeros((y_m.shape[0],2), dtype = 'int32')
        median_pts [:,0] = y_m; median_pts[:,1] = x_m; 
        median_pts *= downscale
    else:    
        print('default')
        #all_pos_outline =  np.column_stack((outline['spatial_1'], outline['spatial_2']))
        p_a, n_a = find_closest_point(pos_a, all_pos_outline, False)
        p_b, n_b = find_closest_point(pos_b, all_pos_outline, False)
        print(p_a)
        print(p_b)
        first_m_point = find_middle_point(p_a, p_b)
        #va = np.array([0,0]); vb = np.array([0,0]); 
        if n_b>n_a:
            all_pos3 = np.delete(all_pos_outline, n_b, 0)
            all_pos3 = np.delete(all_pos3, n_a, 0)
        else:
            all_pos3 = np.delete(all_pos_outline, n_a, 0)
            all_pos3 = np.delete(all_pos3, n_b, 0)
        median_pts, point_a,vector_a,point_b,vector_b, all_pos_crop = iterate_through_points(p_a,v_a, p_b, v_b, all_pos3, first_m_point, shorten_median_prc)
    return median_pts
    
def parabola(x, x0, A, B, y0):
    return A*(x-x0)**2 + B*(x-x0) + y0

def spline(x, x0, A, B, y0):
    return A*(x-x0) + B*(x-x0)**0.5 + y0

def straighten_curve_fit(median_pts, roi_pts, xlim, out_folder, save_images, name, fit_curve=True):
    #print('xlim')
    #print(xlim)
    dx = xlim[1] - xlim[0]
    list_of_x = np.arange(xlim[0]-int(dx/2), xlim[1]+int(dx/2)) #just to make sure we fit slightly out of the xlim
    if fit_curve:
        pop, pcov = curve_fit(parabola, median_pts[:,0], median_pts[:,1])
        #pop, pcov = curve_fit(spline, median_pts[:,0], median_pts[:,1])
        fit = parabola(list_of_x, *pop)
        x_c = find_parabola_xc(pop[0],pop[1],pop[2],pop[3])
        print(x_c)
        n1 = np.where(np.round(list_of_x)==int(x_c)); n1 = int(n1[0])
        #fit = spline(list_of_x, *pop)
    else:
        list_of_x = median_pts[:,0]; fit = median_pts[:,1]
        x_c = np.array([13000, 34000])
        p, n1 = find_closest_point(x_c, median_pts, True)
        
    if save_images:
        plot_and_save_fit(list_of_x, fit, median_pts, out_folder,roi_pts, name)
    
    
    #roi_pts = np.column_stack((ROI['spatial_1'], ROI['spatial_2']))
    str_pts_set = straighten_curve(roi_pts[:,0], roi_pts[:,1], list_of_x, fit, n1)
    return str_pts_set

def plot_and_save(ROI_part_pts, median_pts, median_str_pts, ROI_str_pts, out_folder, name):
    #ROI_part_pts = np.column_stack((ROI_part['spatial_1'], ROI_part['spatial_2']))
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    axs[0].scatter(ROI_part_pts[:,0], ROI_part_pts[:,1], s= 4, color = 'blue', alpha = 0.5)
    axs[0].scatter(median_pts[:,0], median_pts[:,1], s= 7, color = 'red')
    plt.axis('equal')

    axs[1].scatter(ROI_str_pts[:,0], ROI_str_pts[:,1], s= 4, color = 'blue')
    axs[1].scatter(median_str_pts[:,0], median_str_pts[:,1], s= 7, color = 'red', alpha = 0.5)
    path_fig = out_folder + '/' + name + '_curve_plot.png'
    plt.axis('equal')

    plt.savefig(path_fig, format="jpg", dpi = 150)

def plot_and_save_fit(list_of_x, fit, median_points, out_folder, roi_pts, name):
    plt.plot(list_of_x, fit, color = 'r')
    plt.scatter(median_points[:,0], median_points[:,1], s= 5, color = 'blue')
    plt.scatter(roi_pts[:,0], roi_pts[:,1], s= 2, color = 'green')
    plt.xlim((np.min(roi_pts[:,0])-1000, np.max(roi_pts[:,0])+1000))
    plt.ylim((np.min(roi_pts[:,1])-1000, np.max(roi_pts[:,1])+1000))
    path_fig2 = out_folder + '/' + name + '_fit_plot.png'
    plt.savefig(path_fig2, format="jpg", dpi = 150)

def rotate_custom_angle(set_pts, point_rotation, angle):
    #set_pts[:,0] -= point_rotation[0]; set_pts[:,1] -= point_rotation[1]; 
    set_pts_rotated = set_pts.copy(); angle *= (np.pi/180)
    set_pts_rotated[:,0] = (set_pts_rotated[:,0]-point_rotation[0])*np.cos(angle) - (set_pts_rotated[:,1]-point_rotation[1])*np.sin(angle) + point_rotation[0]
    set_pts_rotated[:,1] = (set_pts_rotated[:,1]-point_rotation[1])*np.cos(angle) + (set_pts_rotated[:,0]-point_rotation[0])*np.sin(angle) + point_rotation[1]
    return set_pts_rotated

def rotate_custom_angle2(point, origin, degrees):
    radians = np.deg2rad(degrees)
    x = point[:,0]; y = point[:,1]
    offset_x, offset_y = origin
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = np.cos(radians)
    sin_rad = np.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
    rot_pts = point.copy(); rot_pts[:,0] = qx; rot_pts[:,1] = qy
    return rot_pts

def merge_roi_median_save(ROI_part, ROI_str_pts, median_str_pts, name, out_folder):
    ROI_part['spatial_1'] = ROI_str_pts[:,0]; ROI_part['spatial_2'] = ROI_str_pts[:,1] 
    '''
    blanc_df = ROI_part.iloc[0:2].copy()
    blanc_df['array_row'] = [np.nan, np.nan]; blanc_df['array_col'] = [np.nan, np.nan];  blanc_df['Annotations'] = ['Median', 'Median'] 
    blanc_df = blanc_df.rename(index = {ROI_part.index[0]:'lol', ROI_part.index[1]:'lol'})
    median_df = blanc_df
    median_df['spatial_1'] = median_str_pts[0:2,0]; median_df['spatial_2'] = median_str_pts[0:2,1]; 
    #median_df['spatial_1'][0] = median_str_pts[0,0]; median_df['spatial_2'][0] = median_str_pts[0,1]
    #median_df['spatial_1'][1] = median_str_pts[1,0]; median_df['spatial_2'][1] = median_str_pts[1,1]
    for i in range(2, median_str_pts.shape[0]):
        row = blanc_df.iloc[0].copy(); row['spatial_1'] = median_str_pts[i,0]; row['spatial_2'] = median_str_pts[i,1]
        median_df = median_df.append(row,ignore_index=False) 
    
    ROI_part_median = pd.concat([ROI_part, median_df]) dont need median
    '''
    ROI_part = ROI_part.reset_index().rename(columns={"index":"label_id"})
    path_csv = out_folder + '/' + name + '_straight_pos.csv'
    ROI_part.to_csv(path_csv)
    
def main(ConfFilePath):
    print('preparing data')
    path_adata, name, x_lim, y_lim, pos_a, pos_b, out_folder, rotate_90, U_shape, save_images, shorten_median_prc, subcrop_x, subcrop_y, downscale, mode, rotate_angle, rotate_c_point, v_a, v_b = ReadConfFile(ConfFilePath)
    #print((np.array(subcrop_x)).ndim)
    #print(subcrop_x[1])
    ROI, outline  = get_roi_outline(path_adata)
    ROI_part, outline_part = crop_roi_outline(ROI, outline, x_lim, y_lim)
    #print(outline_part.shape)
    if not subcrop_x == 'None':
        print('CROP')
        if (np.array(subcrop_x)).ndim == 1:
            ROI_part, outline_part = subcrop_roi_outline(ROI_part, outline_part, subcrop_x, subcrop_y)
        else:
            for n_crop_region in range((np.array(subcrop_x)).shape[0]):
                ROI_part, outline_part = subcrop_roi_outline(ROI_part, outline_part, subcrop_x[n_crop_region], subcrop_y[n_crop_region])
        
    
    #if needed - change x to y and y to x
    roi_pts = np.column_stack((ROI_part['spatial_1'], ROI_part['spatial_2']))
    outline_part_pts = np.column_stack((outline_part['spatial_1'], outline_part['spatial_2']))
    if rotate_90: 
        roi_pts[:, 0], roi_pts[:, 1] = roi_pts[:, 1], roi_pts[:, 0].copy()
        outline_part_pts[:, 0], outline_part_pts[:, 1] = outline_part_pts[:, 1], outline_part_pts[:, 0].copy()
        pos_a[0], pos_a[1] = pos_a[1], pos_a[0].copy()
        pos_b[0], pos_b[1] = pos_b[1], pos_b[0].copy()
        v_a[0], v_a[1] = v_a[1], v_a[0].copy()
        v_b[0], v_b[1] = v_b[1], v_b[0].copy()
    print('getting median curve')
    median_pts = define_median_modes(outline_part_pts, pos_a, pos_b, v_a, v_b, shorten_median_prc, out_folder, name, downscale, mode)
    
    if rotate_angle!=0:
        roi_pts = rotate_custom_angle2(roi_pts, np.array(rotate_c_point), rotate_angle)
        median_pts = rotate_custom_angle2(median_pts, np.array(rotate_c_point), rotate_angle)
        x_lim = np.array([np.min(roi_pts[:,0]), np.max(roi_pts[:,0])])
        #x_lim = x_lim - np.array(rotate_c_point[0])
        #y_lim = y_lim - np.array(rotate_c_point[1])

    
    print('straighten ROI points')
    if rotate_90:
        ROI_str_pts = straighten_curve_fit(median_pts, roi_pts, y_lim, out_folder, save_images, name)
        median_str_pts = straighten_curve_fit(median_pts, median_pts, y_lim, out_folder, save_images, name)
    else:
        ROI_str_pts = straighten_curve_fit(median_pts, roi_pts, x_lim, out_folder, save_images, name)
        median_str_pts = straighten_curve_fit(median_pts, median_pts, x_lim, out_folder, save_images, name)
    
                  
    if save_images:
        plot_and_save(roi_pts, median_pts, median_str_pts, ROI_str_pts, out_folder, name)           
                      
    #saving new spot positions and median
    merge_roi_median_save(ROI_part, ROI_str_pts, median_str_pts, name, out_folder)
    
   
    
if __name__ == "__main__":
    fire.Fire(main)     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
