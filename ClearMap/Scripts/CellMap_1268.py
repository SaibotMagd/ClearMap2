#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CellMap
=======

This script is the main pipeline to analyze immediate early gene expression 
data from iDISCO+ cleared tissue [Renier2016]_.

See the :ref:`CellMap tutorial </CellMap.ipynb>` for a tutorial and usage.


.. image:: ../Static/cell_abstract_2016.jpg
   :target: https://doi.org/10.1016/j.cell.2020.01.028
   :width: 300

.. figure:: ../Static/CellMap_pipeline.png

  iDISCO+ and ClearMap: A Pipeline for Cell Detection, Registration, and 
  Mapping in Intact Samples Using Light Sheet Microscopy.


References
----------
.. [Renier2016] `Mapping of brain activity by automated volume analysis of immediate early genes. Renier* N, Adams* EL, Kirst* C, Wu* Z, et al. Cell. 2016 165(7):1789-802 <https://doi.org/10.1016/j.cell.2016.05.007>`_
"""
__author__    = 'Christoph Kirst <christoph.kirst.ck@gmail.com>'
__license__   = 'GPLv3 - GNU General Pulic License v3 (see LICENSE)'
__copyright__ = 'Copyright © 2020 by Christoph Kirst'
__webpage__   = 'http://idisco.info'
__download__  = 'http://www.github.com/ChristophKirst/ClearMap2'

if __name__ == "__main__":
     
  #%%############################################################################
  ### Initialization 
  ###############################################################################
  import time
  starttime = time.time()
  import numpy as np
  
  #ClearMap path
  import sys
  sys.path.append('/home/user/anaconda3/clearmap/ClearMap2')
  #%% Initialize workspace
  
  import pyqtgraph as pg
  print(pg.QtCore.PYQT_VERSION_STR)
  
  from ClearMap.Environment import *  #analysis:ignore
  
  #directories and files
  directory = '/home/user/anaconda3/clearmap/ClearMap2/ClearMap/Tests/Data/CellMap_Example/haloperidol/1268'
  
  expression_raw      = '150819_0_8X-cfos_14-22-33/14-22-33_0_8X-cfos_UltraII_C00_xyz-Table Z<Z,4>.ome.npy'           
  expression_auto     = '150819_0_8X-autofluo_15-30-26/15-30-26_0_8X-autofluo_UltraII_C00_xyz-Table Z<Z,4>.ome.npy'  
  
  ws = wsp.Workspace('CellMap', directory=directory);
  ws.update(raw=expression_raw, autofluorescence=expression_auto)
  ws.info()
  
  ws.debug = False
  
  resources_directory = settings.resources_path
  
  #p3d.plot(ws.filename('raw'))
  #p3d.plot(ws.file_list('raw')[0:2])
  #p3d.plot([ws.file_list('raw')[0:2]])
  """
  """
  #%% Initialize alignment 
  
  #init atlas and reference files
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(0,228)), orientation=(1,2,3),
      overwrite=False, verbose=True);
  
  #alignment parameter files    
  align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bspline.txt')
  
  
  #%%############################################################################
  ### Data conversion
  ############################################################################### 
  
  #%% Convert raw data to npy file     
  """               
  source = ws.source('raw');
  sink   = ws.filename('stitched')
  io.delete_file(sink)
  io.convert(source, sink, processes=None, verbose=True);
  """
  io.mhd.write_header_from_source(ws.filename('stitched'))
  #%%############################################################################
  ### Resampling and atlas alignment 
  ###############################################################################
        
  #%% Resample 
             
  resample_parameter = {
      "source_resolution" : (4.0625, 4.0625, 3),
      "sink_resolution"   : (25,25,25),
      "processes" : 4,
      "verbose" : True,             
      };
  
  io.delete_file(ws.filename('resampled'))
  
  res.resample(ws.filename('stitched'), sink=ws.filename('resampled'), **resample_parameter)
  
  #%%
  #p3d.plot(ws.filename('resampled'))
  
  #%% Resample autofluorescence
      
  resample_parameter_auto = {
      "source_resolution" : (4.0625, 4.0625, 3),
      "sink_resolution"   : (25,25,25),
      "processes" : 4,
      "verbose" : True,                
      };    
  
  res.resample(ws.filename('autofluorescence'), sink=ws.filename('resampled', postfix='autofluorescence'), **resample_parameter_auto)
  
  #p3d.plot([ws.filename('resampled'), ws.filename('resampled', postfix='autofluorescence')])
  #%% Aignment - resampled to autofluorescence (2016er testset 4.55min)
  
  # align the two channels
  align_channels_parameter = {            
      #moving and reference images
      "moving_image" : ws.filename('resampled', postfix='autofluorescence'),
      "fixed_image"  : ws.filename('resampled'),
      
      #elastix parameter files for alignment
      "affine_parameter_file"  : align_channels_affine_file,
      "bspline_parameter_file" : None,
      
      #directory of the alig'/home/nicolas.renier/Documents/ClearMap_Ressources/Par0000affine.txt',nment result
      "result_directory" :  ws.filename('resampled_to_auto')
      }; 
  
  elx.align(**align_channels_parameter);
  
  #%% Alignment - autoflourescence to reference
  
  # align autofluorescence to reference
  align_reference_parameter = {            
      #moving and reference images
      "moving_image" : reference_file,
      "fixed_image"  : ws.filename('resampled', postfix='autofluorescence'),
                                   
      
      #elastix parameter files for alignment
      "affine_parameter_file"  :  align_reference_affine_file,
      "bspline_parameter_file" :  align_reference_bspline_file,
      #directory of the alignment result
      "result_directory" :  ws.filename('auto_to_reference')
      };
  
  elx.align(**align_reference_parameter);
  
  #%%############################################################################
  ### Create test data
  ###############################################################################
  
  #%% Crop test data 
  #print some descriptives about the datasets
  s = ws.source('raw')
  print(s)
  print(s.shape)

  s = ws.source('autofluorescence')
  print(s)
  print(s.shape)
  """
  #select sublice for testing the pipeline
  slicing = (slice(1500,1800),slice(1500,1800),slice(1000,1200));
  ws.create_debug('stitched', slicing=slicing);
  ws.debug = True; 
  
  #p3d.plot(ws.filename('stitched'))
  """    
  #p3d.plot(ws.filename('autofluorescence'))
  #%%############################################################################
  ### Cell detection
  ###############################################################################
 
  #%% Cell detection:
  
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter['illumination'] = None;
  cell_detection_parameter['background'] = None;
  cell_detection_parameter['intensity_detection']['measure'] = ['source'];
  cell_detection_parameter['shape_detection']['threshold'] = 700;
  
  io.delete_file(ws.filename('cells', postfix='maxima'))
  cell_detection_parameter['maxima_detection']['save'] = ws.filename('cells', postfix='maxima')
  
  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
      processes = 6, # 'serial',
      size_max = 20, #100, #35,
      size_min = 11,# 30, #30,
      overlap  = 10, #32, #10,
      verbose = True
      )
  
  cells.detect_cells(ws.filename('stitched'), ws.filename('cells', postfix='raw'),
                     cell_detection_parameter=cell_detection_parameter, 
                     processing_parameter=processing_parameter)
  
  #%% visualization
  #p3d.plot([[ws.filename('stitched'), ws.filename('cells', postfix='maxima')]])
  
  #%%
  #coordinates = np.hstack([ws.source('cells', postfix='raw')[c][:,None] for c in 'xyz']);
  #p = p3d.list_plot_3d(coordinates)
  #p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))
    
  
  #%% Cell statistics
  """
  source = ws.source('cells', postfix='raw')
  
  plt.figure(1); plt.clf();
  names = source.dtype.names;
  nx,ny = p3d.subplot_tiling(len(names));
  for i, name in enumerate(names):
    plt.subplot(nx, ny, i+1)
    plt.hist(source[name]);
    plt.title(name)
  plt.tight_layout();
  print(name)
  print(len(source[name]))
  """
  #%% Filter cells (in 4um resolution 20, in hres. increase it to 60-80)
  
  thresholds = {
      'source' : None,
      'size'   : (20,None)
      }

  cells.filter_cells(source = ws.filename('cells', postfix='raw'), 
                     sink = ws.filename('cells', postfix='filtered'), 
                     thresholds=thresholds);
  
  #import ClearMap.Utils.HierarchicalDict as hdict
  #hdict.pprint(cells.default_cell_detection_parameter)  
  #%% Visualize

  #p3d.plot([source, sink])
  
  
  coordinates = np.array([ws.source('cells', postfix='filtered')[c] for c in 'xyz']).T;
  p = p3d.list_plot_3d(coordinates, color=(1,0,0,0.5), size=10)
  #p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))


  #%%############################################################################
  ### Cell atlas alignment and annotation
  ###############################################################################
  
  #%% Cell alignment
  source = ws.source('cells', postfix='filtered')
  
  def transformation(coordinates):
    coordinates = res.resample_points(
                    coordinates, sink=None, orientation=None, 
                    source_shape=io.shape(ws.filename('stitched')), 
                    sink_shape=io.shape(ws.filename('resampled')));
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('resampled_to_auto'), 
                    binary=True, indices=False);
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('auto_to_reference'),
                    binary=True, indices=False);
        
    return coordinates;
    
  
  coordinates = np.array([source[c] for c in 'xyz']).T;
  
  coordinates_transformed = transformation(coordinates);
  
  #%% Cell annotation  
  label = ano.label_points(coordinates_transformed, key='order');
  names = ano.convert_label(label, key='order', value='name');
  
  #%% Save results
  
  coordinates_transformed.dtype=[(t,float) for t in ('xt','yt','zt')]
  label = np.array(label, dtype=[('order', int)]);
  names = np.array(names, dtype=[('name', 'a256')])
  
  #print(source[0:10])
  print(coordinates_transformed[0:10])
  import numpy.lib.recfunctions as rfn
  cells_data = rfn.merge_arrays([source[:], coordinates_transformed, label, names], flatten=True, usemask=False)
  
  io.write(ws.filename('cells'), cells_data)
  
  
  #%%############################################################################
  ### Cell csv generation for external analysis
  ###############################################################################
  
  #%% CSV export  
  source = ws.source('cells');
  header = ':'.join([h[:] for h in source.dtype.names]);
  fmt = '%.18e:%.18e:%.18e:%.18e:%.18e:%.18e:%.18e:%.18e:%.18e:%s';
  np.savetxt(ws.filename('cells', extension='csv'), source[:], fmt=fmt, header=header, delimiter=':', comments='')  
  
    
  #%% ClearMap 1.0 export
  """
  source = ws.source('cells');
  
  clearmap1_format = {'points' : ['x', 'y', 'z'], 
                      'points_transformed' : ['xt', 'yt', 'zt'],
                      'intensities' : ['source', 'dog', 'background', 'size']}
  
  for filename, names in clearmap1_format.items():
    sink = ws.filename('cells', postfix=['ClearMap1', filename]);
    data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan) for name in names]);
    io.write(sink, data);
  """
  #%%############################################################################
  ### Voxelization - cell density
  ###############################################################################
  
  source = ws.source('cells')
  
  coordinates = np.array([source[n] for n in ['xt','yt','zt']]).T;
  intensities = source['source'];
  
  #%% Unweighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file), 
        dtype = None, 
        weights = None,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='counts'), **voxelization_parameter);
  
  
  #%% 
  
  #p3d.plot(ws.filename('density', postfix='counts'))
  
  
  #%% Weighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file),
        dtype = None, 
        weights = intensities,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='intensities'), **voxelization_parameter);
  
  #%%
  
  #p3d.plot([ws.filename('density', postfix='intensities'),ws.filename('density', postfix='counts')])
  runtime = time.time()-starttime
  print(runtime)