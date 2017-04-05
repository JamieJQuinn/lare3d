#!/usr/bin/env python3

import sys

DeleteAllPlots()

OpenDatabase("localhost:/maths/scratch/lare3d_runs/MT_switch_400/0005.sdf", 0)
AddPlot("Pseudocolor", "Fluid/Anisotropic_Viscous_Heating", 1, 0)
AddOperator("Slice", 0)

DrawPlots()

AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes2D.visible = 1
AnnotationAtts.axes2D.autoSetTicks = 0
AnnotationAtts.axes2D.xAxis.title.userUnits = 1
AnnotationAtts.axes2D.xAxis.title.units = ""
AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = -3
AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 3
AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.1
AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 1
AnnotationAtts.axes2D.yAxis.title.userUnits = 1
AnnotationAtts.axes2D.yAxis.title.units = ""
AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = -3
AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 3
AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.1
AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 1
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
SetAnnotationAttributes(AnnotationAtts)

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 1e-06
PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "Spectral"
PseudocolorAtts.invertColorTable = 1
SetPlotOptions(PseudocolorAtts)

legend = GetAnnotationObject(GetPlotList().GetPlots(0).plotName)
legend.SetDrawMinMax(0)
legend.SetDrawTitle(0)
legend.SetNumTicks(3)
legend.SetNumberFormat('% -9.4g')

SaveWindow()

sys.exit()
