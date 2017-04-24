#!/usr/bin/env python3

import sys
import re

def save_fig(data_fname, fig_fname, is_iso, visc_mode, for_poster=False):
  OpenDatabase("localhost:"+data_fname)

  DeleteAllPlots()

  if visc_mode:
    legend_max = 1e-6
  else:
    legend_max = 1e-5

  if is_iso:
    AddPlot("Pseudocolor", "Fluid/Isotropic_Viscous_Heating", 1, 0)
    fig_fname += "_iso"
  else:
    AddPlot("Pseudocolor", "Fluid/Anisotropic_Viscous_Heating", 1, 0)
    fig_fname += "_aniso"

  AddOperator("Slice", 0)

  DrawPlots()

  AnnotationAtts = AnnotationAttributes()
  if for_poster:
    AnnotationAtts.axes2D.visible = 0
  else:
    AnnotationAtts.axes2D.visible = 1
  AnnotationAtts.axes2D.autoSetTicks = 0
  AnnotationAtts.axes2D.xAxis.title.userUnits = 1
  AnnotationAtts.axes2D.xAxis.title.units = ""
  AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = -3
  AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 3
  AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.1
  AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 3
  AnnotationAtts.axes2D.yAxis.title.userUnits = 1
  AnnotationAtts.axes2D.yAxis.title.units = ""
  AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = -3
  AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 3
  AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.1
  AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 3
  AnnotationAtts.userInfoFlag = 0
  AnnotationAtts.databaseInfoFlag = 0
  AnnotationAtts.timeInfoFlag = 1
  if for_poster:
    AnnotationAtts.legendInfoFlag = 0
  else:
    AnnotationAtts.legendInfoFlag = 1
  SetAnnotationAttributes(AnnotationAtts)

  PseudocolorAtts = PseudocolorAttributes()
  PseudocolorAtts.minFlag = 1
  PseudocolorAtts.min = 0
  PseudocolorAtts.maxFlag = 1
  PseudocolorAtts.max = legend_max
  PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
  PseudocolorAtts.colorTableName = "Spectral"
  PseudocolorAtts.invertColorTable = 1
  SetPlotOptions(PseudocolorAtts)

  legend = GetAnnotationObject(GetPlotList().GetPlots(0).plotName)
  legend.SetDrawMinMax(0)
  legend.SetDrawTitle(0)
  legend.SetNumTicks(3)
  legend.SetNumberFormat('% -9.4g')

  swatts = SaveWindowAttributes()
  swatts.family = 0
  swatts.format = swatts.PNG
  swatts.width = 2048
  swatts.height = 2048
  swatts.fileName = fig_fname
  SetSaveWindowAttributes(swatts)

  if for_poster:
    viewatts = View2DAttributes()
    viewatts.windowCoords = (-3, 3, -3, 3)
    viewatts.viewportCoords = (-3, 3, -3, 3)
    SetView2D(viewatts)

  SaveWindow()

# Get files to generate figs from
with open("filename_list", 'r') as fp:
  filenames = fp.read().split()

for filename in filenames:
  # scrape info from filenames
  n_grid_points, visc_mode = re.search('/(\d+)(-[bs])?/', filename).groups()
  t_step = re.search('/(\d+).sdf', filename).groups()[0]

  # Form output filenames
  if visc_mode:
    fig_fname = n_grid_points + visc_mode + "_" + t_step
  else:
    fig_fname = n_grid_points + "_" + t_step

  # Is this for a poster?
  for_poster = True

  # Save isotropic
  save_fig(filename, fig_fname, True, visc_mode, for_poster)
  if visc_mode:
    # Save anisotropic
    save_fig(filename, fig_fname, False, visc_mode, for_poster)

sys.exit()
