#!/usr/bin/python
#
# Utility for visualizing quality of peak finding.
#
# Hazen 05/13
#

import numpy
import os
import sys
from PyQt4 import QtCore, QtGui

import sa_library.datareader as datareader
import sa_library.readinsight3 as readinsight3

import qtRangeSlider

# UIs.
import visualizer_ui as visualizerUi


#
# Handle Info Table.
#
class InfoTable(QtGui.QWidget):

    def __init__(self, table_widget, specs, parent = None):
        QtGui.QWidget.__init__(self, parent)

        self.specs = specs
        self.table_widget = table_widget

        # setup info table
        self.table_widget.setRowCount(len(specs))
        self.table_widget.setColumnCount(2)

        for i, spec in enumerate(specs):
            widget = QtGui.QTableWidgetItem(spec[0])
            widget.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_widget.setItem(i,0,widget)
            widget = QtGui.QTableWidgetItem("na")
            widget.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_widget.setItem(i,1,widget)

        self.table_widget.resizeColumnToContents(0)

    def update(self, mol_values):
        if mol_values:
            for i, spec in enumerate(self.specs):
                val = "na"
                if (spec[1] in mol_values):
                    val = str(mol_values[spec[1]])
                self.table_widget.item(i,1).setText(val)


#
# Molecule item for the graphics scene.
#
class MoleculeItem(QtGui.QGraphicsEllipseItem):

    def __init__(self, x, y, w, h, type):

        if (type == "i3"):
            self.pen = QtGui.QPen(QtGui.QColor(255,0,0))
            self.pen.setWidthF(0.2)
            self.sel_pen = QtGui.QPen(QtGui.QColor(255,100,200))
            self.sel_pen.setWidthF(0.2)
        else:
            self.pen = QtGui.QPen(QtGui.QColor(0,255,0))
            self.pen.setWidthF(0.6)
            self.sel_pen = QtGui.QPen(QtGui.QColor(100,255,200))
            self.sel_pen.setWidthF(0.6)

        x = x - 0.5*w - 0.5
        y = y - 0.5*h - 0.5
        QtGui.QGraphicsEllipseItem.__init__(self, x, y, w, h)
        self.setPen(self.pen)

    def setMarked(self, marked):
        if marked:
            self.setPen(self.sel_pen)
        else:
            self.setPen(self.pen)


#
# Handle molecule list.
#
class MoleculeList():

    def __init__(self, filename, frame_sx, frame_sy, type):

        self.i3_bin = readinsight3.I3Reader(filename)

        self.frame_sx = frame_sx
        self.frame_sy = frame_sy
        self.last_frame = -1
        self.last_i = 0
        self.mol_items = []
        self.type = type

    def adjustForSetup(self, nm_per_pixel):
        ax = self.data['ax']
        w = self.data['w']
        self.wx = numpy.sqrt(w*w/ax)/nm_per_pixel
        self.wy = numpy.sqrt(w*w*ax)/nm_per_pixel
        #self.x = self.frame_sx - self.data['y'] + 0.5
        #self.y = self.frame_sy - self.data['x'] + 0.5
        self.x = self.data['x']
        self.y = self.data['y']

    def createMolItems(self, frame_number, nm_per_pixel):

        # Only load new localizations if this is a different frame.
        if (frame_number != self.last_frame):
            self.last_frame = frame_number
            self.data = self.i3_bin.getMoleculesInFrameRange(frame_number, frame_number+1)
            self.last_i = 0

        self.adjustForSetup(nm_per_pixel)
        self.mol_items = []
        for i in range(self.x.size):
            self.mol_items.append(MoleculeItem(self.x[i],
                                               self.y[i],
                                               3.0*self.wx[i],
                                               3.0*self.wy[i],
                                               self.type))
        return self.mol_items

    def getClosest(self, px, py):

        if (self.x.size > 0):

            # find the one nearest to px, py.
            dx = self.x - px - 0.5
            dy = self.y - py - 0.5
            dist = dx*dx+dy*dy
            i = numpy.argmin(dist)

            # unmark old item, mark new item
            self.mol_items[self.last_i].setMarked(False)
            self.mol_items[i].setMarked(True)
            self.last_i = i

            # create a dictionary containing the data for this molecule.
            vals = {}
            for field in self.data.dtype.names:
                temp = self.data[field]
                vals[field] = temp[i]

            # add some other useful properties.
            vals['x'] = self.x[i]
            vals['y'] = self.y[i]
            vals['wx'] = self.wx[i]
            vals['wy'] = self.wy[i]

            return vals

        else:
            return False

    def getField(self, field):
        field = str(field)
        if field in self.data.dtype.names:
            return self.data[field]
        elif (field == "wx"):
            return self.wx
        elif (field == "wy"):
            return self.wy
        else:
            return numpy.array([])


#
# Movie view window.
#
class MovieView(QtGui.QGraphicsView):

    #key_press = QtCore.pyqtSignal(object)
    mouse_press = QtCore.pyqtSignal(float, float, name='mousePress')

    def __init__(self, parent, xyi_label):
        QtGui.QGraphicsView.__init__(self, parent)

        # Class variables.
        self.data = False
        self.image = False
        self.margin = 128.0
        self.xyi_label = xyi_label
        self.zoom_in = 1.2
        self.zoom_out = 1.0/self.zoom_in

        # UI initializiation.
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QtCore.QSize(200, 200))

        # Scene initialization.
        self.scene = QtGui.QGraphicsScene()
        self.setScene(self.scene)
        self.setMouseTracking(True)
        self.setRenderHint(QtGui.QPainter.Antialiasing + QtGui.QPainter.SmoothPixmapTransform)
        self.scale(2.0, 2.0)

        # Tooltip.
        self.setToolTip("Advance frame +- 1 <,><.>\nAdvance frame +-200<k><l>\n<Home>first frame\n<End>last frame")

    def keyPressEvent(self, event):
        # This allows keyboard scrolling to work.
        QtGui.QGraphicsView.keyPressEvent(self, event)

        # This allows us to scroll through the movie.
        #self.key_press.emit(event)

    def mouseMoveEvent(self, event):
        pointf = self.mapToScene(event.pos())
        x = pointf.x()
        y = pointf.y()
        i = 0
        if (type(self.data) == type(numpy.array([]))):
            [sx, sy] = self.data.shape
            if ((x>=0) and (x<sx) and (y>=0) and (y<sy)):
                i = int(self.data[int(y), int(x)])
        self.xyi_label.setText("{0:.2f}, {1:.2f}, {2:d}".format(x + 0.5, y + 0.5, i))

    def mousePressEvent(self, event):
        if (event.button() == QtCore.Qt.LeftButton):
            pointf = self.mapToScene(event.pos())
            self.mouse_press.emit(pointf.x(), pointf.y())

    def newFrame(self, frame, multi_molecules, i3_molecules, fmin, fmax):
        self.scene.clear()

        ## process image
        # save image
        self.data = frame.copy()

        # scale image.
        frame = 255.0*(frame-fmin)/(fmax-fmin)
        frame[(frame > 255.0)] = 255.0
        frame[(frame < 0.0)] = 0.0
        
        # convert to QImage.
        frame = numpy.ascontiguousarray(frame.astype(numpy.uint8))
        h, w = frame.shape
        frame_RGB = numpy.zeros((frame.shape[0], frame.shape[1], 4), dtype = numpy.uint8)
        frame_RGB[:,:,0] = frame
        frame_RGB[:,:,1] = frame
        frame_RGB[:,:,2] = frame
        frame_RGB[:,:,3] = 255

        self.image = QtGui.QImage(frame_RGB.data, w, h, QtGui.QImage.Format_RGB32)
        self.image.ndarray1 = frame
        self.image.ndarray2 = frame_RGB
    
        # add to scene
        self.scene.addPixmap(QtGui.QPixmap.fromImage(self.image))

        # add 3D-DAOSTORM localizations
        for loc in multi_molecules:
            self.scene.addItem(loc)

        # add Insight3 localizations
        for loc in i3_molecules:
            self.scene.addItem(loc)

    def wheelEvent(self, event):
        if event.delta() > 0:
            self.zoomIn()
        else:
            self.zoomOut()

    def zoomIn(self):
        self.scale(self.zoom_in, self.zoom_in)

    def zoomOut(self):
        self.scale(self.zoom_out, self.zoom_out)


#
# Main window
#
class Window(QtGui.QMainWindow):
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)

        # variables
        self.cur_frame = 0
        self.directory = ""
        self.film_l = 0
        self.film_x = 255
        self.film_y = 255
        self.i3_list = False
        self.movie_file = False
        self.multi_list = False
        self.settings = QtCore.QSettings("Zhuang Lab", "visualizer")

        self.locs_display_timer = QtCore.QTimer(self)
        self.locs_display_timer.setInterval(100)
        self.locs_display_timer.setSingleShot(True)
        self.locs_display_timer.timeout.connect(self.handleLocsDisplayTimer)

        # ui setup
        self.ui = visualizerUi.Ui_MainWindow()
        self.ui.setupUi(self)

        # initialize info tables
        self.multi_table = InfoTable(self.ui.multiTableWidget,
                                     [["x", "x", "float"],
                                      ["y", "y", "float"],
                                      ["z", "z", "float"],
                                      ["height", "h", "float"],
                                      ["width-x", "wx", "float"],
                                      ["width-y", "wy", "float"],
                                      ["background", "bg", "float"],
                                      ["sum", "a", "float"],
                                      ["fit error", "i", "float"],
                                      ["status", "fi", "int"]])

        self.i3_table = InfoTable(self.ui.i3TableWidget,
                                  [["x", "x", "float"],
                                   ["y", "y", "float"],
                                   ["z", "z", "float"],
                                   ["height", "h", "float"],
                                   ["area (fit)", "a", "float"],
                                   ["width-x", "wx", "float"],
                                   ["width-y", "wy", "float"],
                                   ["background", "bg", "float"],
                                   ["sum", "i", "float"],
                                   ["category", "c", "int"],
                                   ["iteration", "fi", "int"],
                                   ["track length", "tl", "int"]])
        
        # initialize movie viewing tab.
        self.movie_view = MovieView(self.ui.movieGroupBox, self.ui.xyiLabel)
        movie_layout = QtGui.QGridLayout(self.ui.movieGroupBox)
        movie_layout.addWidget(self.movie_view)
        self.movie_view.show()
        #self.movie_view.key_press.connect(self.keyPressEvent)
        self.movie_view.mouse_press.connect(self.updateInfo)

        # initialize range slider.
        self.rangeSlider = qtRangeSlider.QVRangeSlider([self.ui.minSpinBox.minimum(),
                                                        self.ui.maxSpinBox.maximum(),
                                                        1.0],
                                                       [self.ui.minSpinBox.value(),
                                                        self.ui.maxSpinBox.value()],
                                                       parent = self.ui.rangeSliderWidget)
        layout = QtGui.QGridLayout(self.ui.rangeSliderWidget)
        layout.addWidget(self.rangeSlider)
        self.rangeSlider.setEmitWhileMoving(True)
        self.rangeSlider.rangeChanged.connect(self.handleRangeChange)

        # signals
        self.ui.actionCapture.triggered.connect(self.capture)
        self.ui.actionLoad_3DDAO_Locs.triggered.connect(self.load3DDAOLocalizations)
        self.ui.actionLoad_Insight3_Locs.triggered.connect(self.loadI3Localizations)
        self.ui.actionLoad_Movie.triggered.connect(self.loadMovie)
        self.ui.actionQuit.triggered.connect(self.quit)
        self.ui.maxSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.minSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.nmPerPixelSpinBox.valueChanged.connect(self.handleNmPerPixelSpinBox)
        self.ui.oriCheckBox.stateChanged.connect(self.handleCheckBox)

        # load settings.
        self.directory = str(self.settings.value("directory", "").toString())
        self.move(self.settings.value("position", QtCore.QPoint(100, 100)).toPoint())
        self.resize(self.settings.value("size", self.size()).toSize())

    def capture(self):
        pixmap = QtGui.QPixmap.grabWidget(self.movie_view.viewport())
        pixmap.save("capture.png")
        print "Capture size:", pixmap.width(), pixmap.height()
        
    def cleanUp(self):
        self.settings.setValue("directory", self.directory)
        self.settings.setValue("position", self.pos())
        self.settings.setValue("size", self.size())

    def closeEvent(self, event):
        self.cleanUp()

    def displayFrame(self, update_locs):
        if self.movie_file:

            # Get the current frame.
            frame = self.movie_file.loadAFrame(self.cur_frame).astype(numpy.float)
            if self.ui.oriCheckBox.isChecked():
                frame = numpy.rot90(numpy.rot90(frame))
            else:
                frame = numpy.transpose(frame)


            # Create the 3D-DAOSTORM molecule items.
            nm_per_pixel = self.ui.nmPerPixelSpinBox.value()
            multi_mols = []
            if update_locs and self.multi_list:
                multi_mols = self.multi_list.createMolItems(self.cur_frame+1, nm_per_pixel)

            # Create the Insight3 molecule items.
            i3_mols = []
            if update_locs and self.i3_list:
                i3_mols = self.i3_list.createMolItems(self.cur_frame+1, nm_per_pixel)

            self.movie_view.newFrame(frame,
                                     multi_mols,
                                     i3_mols,
                                     self.ui.minSpinBox.value(),
                                     self.ui.maxSpinBox.value())

    def handleCheckBox(self, value):
        self.displayFrame(True)

    def handleLocsDisplayTimer(self):
        self.displayFrame(True)

    def handleMaxMinSpinBox(self, value):
        self.rangeSlider.setValues([self.ui.minSpinBox.value(),
                                    self.ui.maxSpinBox.value()])

    def handleNmPerPixelSpinBox(self, value):
        self.displayFrame(True)

    def handleRangeChange(self, range_min, range_max):
        self.ui.minSpinBox.setValue(range_min)
        self.ui.maxSpinBox.setValue(range_max)
        self.displayFrame(False)
        self.locs_display_timer.start()

    def incCurFrame(self, amount):
        self.cur_frame += amount
        if (self.cur_frame < 0):
            self.cur_frame = 0
        if (self.cur_frame >= self.film_l):
            self.cur_frame = self.film_l - 1
        if self.movie_file:
            self.ui.frameLabel.setText("frame " + str(self.cur_frame+1) + " (" + str(self.film_l) + ")")
            self.displayFrame(False)
            self.locs_display_timer.start()

    def keyPressEvent(self, event):
        if (event.key() == QtCore.Qt.Key_End):
            self.incCurFrame(self.film_l)
        if (event.key() == QtCore.Qt.Key_Home):
            self.incCurFrame(-self.film_l)
        if (event.key() == QtCore.Qt.Key_Comma):
            self.incCurFrame(-1)
        if (event.key() == QtCore.Qt.Key_Period):
            self.incCurFrame(1)
        if (event.key() == QtCore.Qt.Key_K):
            self.incCurFrame(-200)
        if (event.key() == QtCore.Qt.Key_L):
            self.incCurFrame(200)

    def load3DDAOLocalizations(self):
        list_filename = str(QtGui.QFileDialog.getOpenFileName(self,
                                                              "Load 3D-DAOSTORM Localization List",
                                                              self.directory,
                                                              "*.bin"))
        if list_filename:
            self.directory = os.path.dirname(list_filename)
            self.multi_list = MoleculeList(list_filename, self.film_x, self.film_y, "3d")
            self.incCurFrame(0)

    def loadI3Localizations(self):
        list_filename = str(QtGui.QFileDialog.getOpenFileName(self,
                                                              "Load Insight3 Localization List",
                                                              self.directory,
                                                              "*.bin"))
        if list_filename:
            self.directory = os.path.dirname(list_filename)
            self.i3_list = MoleculeList(list_filename, self.film_x, self.film_y, "i3")
            self.incCurFrame(0)

    def loadMovie(self):
        movie_filename = str(QtGui.QFileDialog.getOpenFileName(self,
                                                               "Load Movie",
                                                               self.directory,
                                                               "*.dax *.spe *.tif"))
        if movie_filename:
            self.directory = os.path.dirname(movie_filename)
            self.movie_file = datareader.inferReader(movie_filename)
            [self.film_x, self.film_y, self.film_l] = self.movie_file.filmSize()
            self.ui.fileLabel.setText(movie_filename)
            self.cur_frame = 0
            self.multi_list = False
            self.incCurFrame(0)

    def quit(self):
        self.close()

    def updateInfo(self, x, y):
        if self.multi_list:
            vals = self.multi_list.getClosest(x, y)
            self.multi_table.update(vals)
        if self.i3_list:
            vals = self.i3_list.getClosest(x, y)
            self.i3_table.update(vals)

    def wheelEvent(self, event):
        if event.delta() > 0:
            self.incCurFrame(1)
        else:
            self.incCurFrame(-1)


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()


#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
