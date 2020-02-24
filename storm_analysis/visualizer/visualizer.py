#!/usr/bin/env python
"""
Utility for visualizing quality of peak finding.

Hazen 05/13
"""

import numpy
import os
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.sa_h5py as saH5Py

import qtRangeSlider

# UIs.
import visualizer_ui as visualizerUi


class InfoTable(QtWidgets.QWidget):
    """
    Handle Info Table.
    """

    def __init__(self, table_widget = None, **kwds):
        super(InfoTable, self).__init__(**kwds)

        self.fields = None
        self.table_widget = table_widget

        # Setup info table.
        self.table_widget.setColumnCount(2)

    def hideFields(self):
        self.table_widget.setRowCount(0)
        
    def showFields(self, fields):
        self.fields = fields
        self.table_widget.setRowCount(len(fields))
        for i, field in enumerate(fields):
            widget = QtWidgets.QTableWidgetItem(field)
            widget.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_widget.setItem(i,0,widget)
            widget = QtWidgets.QTableWidgetItem("na")
            widget.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_widget.setItem(i,1,widget)
        self.table_widget.resizeColumnToContents(0)

    def update(self, properties):
        if bool(properties):
            for i, field in enumerate(self.fields):
                val = "na"
                if (field in properties):
                    if isinstance(properties[field], numpy.int32):
                        val = str(properties[field])
                    else:
                        val = "{0:.3f}".format(properties[field])
                self.table_widget.item(i,1).setText(val)


class MoleculeItem(QtWidgets.QGraphicsEllipseItem):
    """
    Molecule item for the graphics scene.
    """
    def __init__(self, x, y, w, h, mtype):

        if (mtype == "l2"):
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
        QtWidgets.QGraphicsEllipseItem.__init__(self, x, y, w, h)
        self.setPen(self.pen)

    def setMarked(self, marked):
        if marked:
            self.setPen(self.sel_pen)
        else:
            self.setPen(self.pen)


class MoleculeList(object):
    """
    Handle molecule list.
    """
    def __init__(self, mtype = None, **kwds):
        super(MoleculeList, self).__init__(**kwds)

        self.fields = []
        self.last_frame = -1
        self.last_index = 0
        self.locs = {}
        self.mol_items = []
        self.mtype = mtype
        self.reader = None

    def clearLocs(self):
        self.last_index = 0
        self.locs = {}

    def createMolItems(self, frame_number, nm_per_pixel):
        self.locs = self.loadLocalizations(frame_number, nm_per_pixel)        
        self.mol_items = []
        if bool(self.locs):
            if "xsigma" in self.locs:
                xs = self.locs["xsigma"]
                ys = self.locs["ysigma"]
            else:
                xs = 1.5*numpy.ones(self.locs["x"].size)
                ys = 1.5*numpy.ones(self.locs["x"].size)
            for i in range(self.locs["x"].size):
                self.mol_items.append(MoleculeItem(self.locs["x"][i] + 1,
                                                   self.locs["y"][i] + 1,
                                                   xs[i]*6.0,
                                                   ys[i]*6.0,
                                                   self.mtype))
        return self.mol_items

    def getClosest(self, px, py):
        """
        Return a dictionary with information about the closest localization
        to px, py.
        """
        vals = {}
        if bool(self.locs) and (self.locs["x"].size > 0):

            # Find the one nearest to px, py.
            dx = self.locs["x"] - px - 0.5
            dy = self.locs["y"] - py - 0.5
            dist = dx*dx+dy*dy
            closest_index = numpy.argmin(dist)

            # Unmark old item, mark new item
            self.mol_items[self.last_index].setMarked(False)
            self.mol_items[closest_index].setMarked(True)
            self.last_index = closest_index

            # Create a dictionary containing the data for this molecule.
            vals = {}
            for field in self.locs:
                vals[field] = self.locs[field][closest_index]

        return vals

    def getFields(self):
        return self.fields
        

class MoleculeListHDF5(MoleculeList):
    """
    Handle HDF5 molecule list.
    """
    def __init__(self, filename = None, **kwds):
        super(MoleculeListHDF5, self).__init__(**kwds)

        self.fields = ["x",
                       "y",
                       "z",
                       "background",
                       "error",
                       "height",
                       "sum",                       
                       "xsigma",
                       "ysigma",
                       "category",
                       "iterations",                       
                       "significance"]

        self.reader = saH5Py.SAH5Py(filename)

    def cleanUp(self):
        self.reader.close(verbose = False)

    def loadLocalizations(self, frame_number, nm_per_pixel):
        locs = self.reader.getLocalizationsInFrame(frame_number)
        if bool(locs):
            if ("xsigma" in locs) and (not "ysigma" in locs):
                locs["ysigma"] = locs["xsigma"]
        return locs
        

class MoleculeListI3(MoleculeList):
    """
    Handle Insight3 molecule list.
    """
    def __init__(self, filename = None, **kwds):
        super(MoleculeListI3, self).__init__(**kwds)

        self.fields = ["x",
                       "y",
                       "z",
                       "background",
                       "error",
                       "height",
                       "sum",
                       "xsigma",
                       "ysigma"]

        self.reader = readinsight3.I3Reader(filename)

        self.n_locs = self.reader.getNumberMolecules()

    def cleanUp(self):
        self.reader.close()

    def loadLocalizations(self, frame_number, nm_per_pixel):
        if (self.n_locs > 0):
            fnum = frame_number + 1
            i3data = self.reader.getMoleculesInFrame(fnum)
            return i3dtype.convertToSAHDF5(i3data, fnum, nm_per_pixel)
        else:
            return {}
            

class MovieView(QtWidgets.QGraphicsView):
    """
    Movie view window.
    """

    #key_press = QtCore.pyqtSignal(object)
    mouse_press = QtCore.pyqtSignal(float, float, name='mousePress')

    def __init__(self, xyi_label = None, **kwds):
        super(MovieView, self).__init__(**kwds)

        # Class variables.
        self.data = False
        self.image = False
        self.margin = 128.0
        self.xyi_label = xyi_label
        self.zoom_in = 1.2
        self.zoom_out = 1.0/self.zoom_in

        # UI initializiation.
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,
                                           QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QtCore.QSize(200, 200))

        # Scene initialization.
        self.scene = QtWidgets.QGraphicsScene()
        self.setScene(self.scene)
        self.setMouseTracking(True)
        self.setRenderHint(QtGui.QPainter.Antialiasing + QtGui.QPainter.SmoothPixmapTransform)
        self.scale(2.0, 2.0)

        # Tooltip.
        self.setToolTip("Advance frame +- 1 <,><.>\nAdvance frame +-200<k><l>\n<Home>first frame\n<End>last frame")

    def keyPressEvent(self, event):
        # This allows keyboard scrolling to work.
        QtWidgets.QGraphicsView.keyPressEvent(self, event)

        # This allows us to scroll through the movie.
        #self.key_press.emit(event)

    def mouseMoveEvent(self, event):
        pointf = self.mapToScene(event.pos())
        x = pointf.x()
        y = pointf.y()
        i = 0
        
        if not isinstance(self.data, bool):
            [sy, sx] = self.data.shape
            if ((x>=0) and (x<sx) and (y>=0) and (y<sy)):
                i = int(self.data[int(y), int(x)])
        self.xyi_label.setText("{0:.2f}, {1:.2f}, {2:d}".format(x - 0.5, y - 0.5, i))

    def mousePressEvent(self, event):
        if (event.button() == QtCore.Qt.LeftButton):
            pointf = self.mapToScene(event.pos())
            self.mouse_press.emit(pointf.x(), pointf.y())

    def newFrame(self, frame, locs1, locs2, fmin, fmax):
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

        # add localizations from file 1
        for loc in locs1:
            self.scene.addItem(loc)

        # add localizations from file 2
        for loc in locs2:
            self.scene.addItem(loc)

    def wheelEvent(self, event):
        if not event.angleDelta().isNull():
            if (event.angleDelta().y() > 0):
                self.zoomIn()
            else:
                self.zoomOut()

            # This blocks propogation to the main window where the
            # same event will cause the frame number to change.
            event.accept()

    def zoomIn(self):
        self.scale(self.zoom_in, self.zoom_in)

    def zoomOut(self):
        self.scale(self.zoom_out, self.zoom_out)


class Window(QtWidgets.QMainWindow):
    """
    Main window.
    """
    
    def __init__(self, directory = None, **kwds):
        super(Window, self).__init__(**kwds)

        # variables
        self.cur_frame = 0
        self.directory = ""
        self.film_l = 0
        self.film_x = 255
        self.film_y = 255
        self.locs1_list = None
        self.locs2_list = None
        self.movie_file = None
        self.settings = QtCore.QSettings("Zhuang Lab", "visualizer")

        self.locs_display_timer = QtCore.QTimer(self)
        self.locs_display_timer.setInterval(100)
        self.locs_display_timer.setSingleShot(True)
        self.locs_display_timer.timeout.connect(self.handleLocsDisplayTimer)

        # ui setup
        self.ui = visualizerUi.Ui_MainWindow()
        self.ui.setupUi(self)

        # initialize info tables
        self.locs1_table = InfoTable(table_widget = self.ui.locs1TableWidget)
        self.locs2_table = InfoTable(table_widget = self.ui.locs2TableWidget)
        
        # initialize movie viewing tab.
        self.movie_view = MovieView(xyi_label = self.ui.xyiLabel, parent = self.ui.movieGroupBox)
        movie_layout = QtWidgets.QGridLayout(self.ui.movieGroupBox)
        movie_layout.addWidget(self.movie_view)
        self.movie_view.show()
        #self.movie_view.key_press.connect(self.keyPressEvent)
        self.movie_view.mouse_press.connect(self.updateInfo)

        self.ui.maxSpinBox.setValue(int(self.settings.value("maximum", 2000)))
        self.ui.minSpinBox.setValue(int(self.settings.value("minimum", 100)))
        self.ui.nmPerPixelSpinBox.setValue(float(self.settings.value("pixel_size", 160)))
        
        # initialize range slider.
        self.rangeSlider = qtRangeSlider.QVRangeSlider([self.ui.minSpinBox.minimum(),
                                                        self.ui.maxSpinBox.maximum(),
                                                        1.0],
                                                       [self.ui.minSpinBox.value(),
                                                        self.ui.maxSpinBox.value()],
                                                       parent = self.ui.rangeSliderWidget)
        layout = QtWidgets.QGridLayout(self.ui.rangeSliderWidget)
        layout.addWidget(self.rangeSlider)
        self.rangeSlider.setEmitWhileMoving(True)
        self.rangeSlider.rangeChanged.connect(self.handleRangeChange)

        # signals
        self.ui.actionCapture.triggered.connect(self.capture)
        self.ui.actionLoad_Locs1.triggered.connect(self.handleLoadLocs1)
        self.ui.actionLoad_Locs2.triggered.connect(self.handleLoadLocs2)
        self.ui.actionLoad_Movie.triggered.connect(self.handleLoadMovie)
        self.ui.actionQuit.triggered.connect(self.quit)
        self.ui.maxSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.minSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.nmPerPixelSpinBox.valueChanged.connect(self.handleNmPerPixelSpinBox)

        # load settings.
        if directory is None:
            self.directory = str(self.settings.value("directory", ""))
        else:
            self.directory = directory
            
        self.move(self.settings.value("position", self.pos()))
        self.resize(self.settings.value("size", self.size()))

    def capture(self):
        pixmap = self.movie_view.grab()
        pixmap.save("capture.png")
        print("Capture size:", pixmap.width(), pixmap.height())
        
    def cleanUp(self):
        for elt in [self.locs1_list, self.locs2_list]:
            if elt is not None:
                elt.cleanUp()

        self.settings.setValue("directory", self.directory)
        self.settings.setValue("maximum", self.ui.maxSpinBox.value())
        self.settings.setValue("minimum", self.ui.minSpinBox.value())
        self.settings.setValue("pixel_size", self.ui.nmPerPixelSpinBox.value())
        self.settings.setValue("position", self.pos())
        self.settings.setValue("size", self.size())

    def closeEvent(self, event):
        self.cleanUp()

    def displayFrame(self, update_locs):
        if self.movie_file:

            # Get the current frame.
            frame = self.movie_file.loadAFrame(self.cur_frame).astype(numpy.float)

            # Create localization list 1 molecule items.
            nm_per_pixel = self.ui.nmPerPixelSpinBox.value()
            locs1 = []
            if update_locs and (self.locs1_list is not None):
                locs1 = self.locs1_list.createMolItems(self.cur_frame, nm_per_pixel)

            # Create localization list 2 molecule items.
            locs2 = []
            if update_locs and (self.locs2_list is not None):
                locs2 = self.locs2_list.createMolItems(self.cur_frame, nm_per_pixel)

            self.movie_view.newFrame(frame,
                                     locs1,
                                     locs2,
                                     self.ui.minSpinBox.value(),
                                     self.ui.maxSpinBox.value())

    def handleLoadLocs1(self):
        list_filename = QtWidgets.QFileDialog.getOpenFileName(self,
                                                              "Load Localization List 1",
                                                              self.directory,
                                                              "*.bin *.hdf5")[0]
        if list_filename:
            if self.locs1_list is not None:
                self.locs1_list.cleanUp()
            self.directory = os.path.dirname(list_filename)
            if saH5Py.isSAHDF5(list_filename):
                self.locs1_list = MoleculeListHDF5(filename = list_filename,
                                                   mtype = "l1")
            else:
                self.locs1_list = MoleculeListI3(filename = list_filename,
                                                 mtype = "l1")
            self.locs1_table.showFields(self.locs1_list.getFields())
            self.incCurFrame(0)

    def handleLoadLocs2(self):
        list_filename = QtWidgets.QFileDialog.getOpenFileName(self,
                                                              "Load Localization List 2",
                                                              self.directory,
                                                              "*.bin *.hdf5")[0]
        if list_filename:
            if self.locs2_list is not None:
                self.locs2_list.cleanUp()
            self.directory = os.path.dirname(list_filename)
            if saH5Py.isSAHDF5(list_filename):
                self.locs2_list = MoleculeListHDF5(filename = list_filename,
                                                   mtype = "l2")
            else:
                self.locs2_list = MoleculeListI3(filename = list_filename,
                                                 mtype = "l2")
            self.locs2_table.showFields(self.locs2_list.getFields())
            self.incCurFrame(0)

    def handleLoadMovie(self):
        movie_filename = QtWidgets.QFileDialog.getOpenFileName(self,
                                                               "Load Movie",
                                                               self.directory,
                                                               "*.dax *.fits *.spe *.tif")[0]
        if movie_filename:            
            self.directory = os.path.dirname(movie_filename)
            self.movie_file = datareader.inferReader(movie_filename)
            [self.film_x, self.film_y, self.film_l] = self.movie_file.filmSize()
            self.ui.fileLabel.setText(movie_filename)
            self.cur_frame = 0

            # Clear molecule lists.
            for elt in [self.locs1_list, self.locs2_list]:
                if elt is not None:
                    elt.cleanUp()
            self.locs1_list = None
            self.locs2_list = None

            # Hide info displays
            self.locs1_table.hideFields()
            self.locs2_table.hideFields()
            
            # Reset view transform.
            self.movie_view.setTransform(QtGui.QTransform())
            
            self.incCurFrame(0)

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
            if (self.locs1_list is not None):
                self.locs1_list.clearLocs()
            if (self.locs2_list is not None):
                self.locs2_list.clearLocs()             
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

    def quit(self):
        self.close()

    def updateInfo(self, x, y):
        if self.locs1_list is not None:
            properties = self.locs1_list.getClosest(x, y)
            self.locs1_table.update(properties)
        if self.locs2_list:
            properties = self.locs2_list.getClosest(x, y)
            self.locs2_table.update(properties)

    def wheelEvent(self, event):
        if not event.angleDelta().isNull():
            if (event.angleDelta().y() > 0):
                self.incCurFrame(1)
            else:
                self.incCurFrame(-1)
                

if (__name__ == "__main__"):
    app = QtWidgets.QApplication(sys.argv)

    directory = None
    if (len(sys.argv) == 2):
        directory = sys.argv[1]
        
    window = Window(directory = directory)
    
    window.show()
    app.exec_()
    app.deleteLater()
    sys.exit()


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
