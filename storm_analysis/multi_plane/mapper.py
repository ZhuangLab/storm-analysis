#!/usr/bin/env python
"""
Utility for identifying the transform between different
channels in a multi-channel data set.

This assumes that each channel is it's own movie and
that the localizations in the channel have been 
identified with a tool like 3D-DAOSTORM.

This only identifies a first order mapping between the
channels.

Hazen 05/17
"""

import numpy
import os
import pickle
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.visualizer.qtRangeSlider as qtRangeSlider
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.sa_h5py as saH5Py


import storm_analysis.multi_plane.mapper_ui as mapperUi


class Channel(object):
    """
    Wrapper for a STORM movie and it's associated localization list.
    """
    def __init__(self, color = None,
                 locs_name = None,
                 movie_name = None,
                 number = None,
                 **kwds):
        super().__init__(**kwds)
        self.color = color
        self.cur_frame = 0
        self.flip_lr = False
        self.flip_ud = False
        self.fr_height = None
        self.fr_width = None
        self.image = None
        self.locs = None
        self.locs_name = locs_name
        self.movie_name = movie_name
        self.number = number
        self.offset_x = 0
        self.offset_y = 0
        self.pixmap = None

        if (locs_name.endswith(".bin")):
            self.locs_reader = LocalizationReaderI3(locs_name)
        else:
            self.locs_reader = LocalizationReaderH5(locs_name)

        self.movie_fp = datareader.inferReader(movie_name)
        self.movie_len = self.movie_fp.filmSize()[2]

    def addFrame(self, scene, frame_no, fmin, fmax):
        """
        Add a frame to a graphics scene.
        """
        frame = self.movie_fp.loadAFrame(self.cur_frame).astype(numpy.float)

        # Scale image.
        frame = 255.0*(frame-fmin)/(fmax-fmin)
        frame[(frame > 255.0)] = 255.0
        frame[(frame < 0.0)] = 0.0
        
        # Convert to QImage.
        frame = numpy.ascontiguousarray(frame.astype(numpy.uint8))
        self.fr_height, self.fr_width = frame.shape
        frame_RGB = numpy.zeros((frame.shape[0], frame.shape[1], 4), dtype = numpy.uint8)
        frame_RGB[:,:,0] = frame
        frame_RGB[:,:,1] = frame
        frame_RGB[:,:,2] = frame
        frame_RGB[:,:,3] = 255

        self.image = QtGui.QImage(frame_RGB.data,
                                  self.fr_width,
                                  self.fr_height,
                                  QtGui.QImage.Format_RGB32)
        self.image.ndarray1 = frame
        self.image.ndarray2 = frame_RGB
    
        # Add to scene
        self.pixmap = QtGui.QPixmap.fromImage(self.image.mirrored(self.flip_lr, self.flip_ud))
        pixmap_item = QtWidgets.QGraphicsPixmapItem(self.pixmap)
        pixmap_item.setOffset(self.offset_x, self.offset_y)
        scene.addItem(pixmap_item)

    def addLocs(self, scene, mappings):
        [x, y] = self.getLocs()
        
        for i in range(x.size):
            loc_item = LocalizationItem(i, x[i], y[i], 6.0, self.color)
            scene.addItem(loc_item)

        m_name = str(self.number) + "_0_"
        if (m_name + "x") in mappings:
            rx = self.locs["x"]
            ry = self.locs["y"]
            mx = mappings[m_name + "x"]
            my = mappings[m_name + "y"]
            for i in range(rx.size):
                xf = mx[0] + mx[1]*rx[i] + mx[2]*ry[i]
                yf = my[0] + my[1]*rx[i] + my[2]*ry[i]
                mapped_item = MappedItem(xf, yf)
                scene.addItem(mapped_item)

    def changeFrame(self, frame_no):
        if (frame_no < 0):
            self.cur_frame = 0
        elif (frame_no >= self.movie_len):
            self.cur_frame = self.movie_len - 1
        else:
            self.cur_frame = frame_no
        return self.cur_frame

    def findNearest(self, xp, yp, max_d):
        [x, y] = self.getLocs()
        if (x.size > 0):
            dx = xp - x[0]
            dy = yp - y[0]
            min_d = dx*dx + dy*dy
            min_i = 0
            for i in range(x.size):
                dx = xp - x[i]
                dy = yp - y[i]
                dd = dx*dx + dy*dy
                if (dd < min_d):
                    min_d = dd
                    min_i = i
            if (min_d < (max_d * max_d)):
                return [[x[min_i], y[min_i]],
                        [self.locs["x"][min_i], self.locs["y"][min_i]]]

    def flipLR(self):
        self.flip_lr = not self.flip_lr

    def flipUD(self):
        self.flip_ud = not self.flip_ud
        
    def getCurrentFrameNumber(self):
        return self.cur_frame

    def getLocs(self):
        self.locs = self.locs_reader.loadLocalizations(self.cur_frame)
        x = self.locs["x"].copy()
        y = self.locs["y"].copy()

        if self.flip_lr:
            x = self.fr_width - x - 1 + self.offset_x
        else:
            x += self.offset_x
            
        if self.flip_ud:
            y = self.fr_height - y - 1 + self.offset_y
        else:
            y += self.offset_y

        return [x, y]
        
    def getMovieLength(self):
        return self.movie_len

    def offsetX(self, inc):
        self.offset_x += inc

    def offsetY(self, inc):
        self.offset_y += inc


class Group(object):
    """
    A group of localizations that represent the position of the
    same object in different channels.
    """
    def __init__(self, **kwds):
        super().__init__(**kwds)
        self.raw_coords = []
        self.trans_coords = []

    def addGroup(self, scene):
        for i in range(len(self.trans_coords) - 1):
            scene.addItem(GroupItem(self.trans_coords[0], self.trans_coords[i+1]))
        
    def addLocalization(self, trans_coords, raw_coords):
        self.trans_coords.append(trans_coords)
        self.raw_coords.append(raw_coords)

    def getSize(self):
        return len(self.raw_coords)

    def getXY(self, channel_number):
        return self.raw_coords[channel_number]

    def isSame(self, other):
        p1 = self.raw_coords[0]
        p2 = other.raw_coords[0]
        return ((p1[0] == p2[0]) and (p1[1] == p2[1]))

        
class GroupItem(QtWidgets.QGraphicsLineItem):
    """
    Group item for a graphics scene.
    """
    def __init__(self, p1, p2):
        x1 = p1[0]+0.5
        y1 = p1[1]+0.5
        x2 = p2[0]+0.5
        y2 = p2[1]+0.5
        QtWidgets.QGraphicsLineItem.__init__(self, x1, y1, x2, y2)
        
        pen = QtGui.QPen(QtGui.QColor(255,255,255))
        pen.setWidthF(1.0)
        self.setPen(pen)
        

class LocalizationItem(QtWidgets.QGraphicsEllipseItem):
    """
    Localization item for a graphics scene.
    """
    def __init__(self, loc_index, x, y, d, color):
        self.loc_index = loc_index

        self.default_pen = QtGui.QPen(color)        
        self.default_pen.setWidthF(0.6)

        x = x - 0.5*d + 0.5
        y = y - 0.5*d + 0.5
        QtWidgets.QGraphicsEllipseItem.__init__(self, x, y, d, d)
        self.setPen(self.default_pen)


class LocalizationReaderH5(object):
    """
    Localization reader for HDF5 files.
    """
    def __init__(self, filename):
        super(LocalizationReaderH5, self).__init__()
        self.reader = saH5Py.SAH5Py(filename)

    def loadLocalizations(self, frame_number):
        return self.reader.getLocalizationsInFrame(frame_number)
    
        
class LocalizationReaderI3(object):
    """
    Localization reader for Insight3 files.
    """
    def __init__(self, filename):
        super(LocalizationReaderI3, self).__init__()
        self.reader = readinsight3.I3Reader(filename)
        
    def loadLocalizations(self, frame_number, nm_per_pixel = 100.0):
        fnum = frame_number + 1
        i3data = self.reader.getMoleculesInFrame(fnum)
        return i3dtype.convertToSAHDF5(i3data, fnum, nm_per_pixel)


class MappedItem(QtWidgets.QGraphicsItem):
    """
    (Mapped) localization item for a graphics scene.
    """
    def __init__(self, x, y):
        super().__init__()
        self.setPos(x + 0.5, y + 0.5)

    def boundingRect(self):
        return QtCore.QRectF(-2, -2, 5, 5)

    def paint(self, painter, option, widget):
        pen = QtGui.QPen(QtGui.QColor(255, 255, 255))
        pen.setWidthF(0.2)
        painter.setPen(pen)
        painter.drawLine(-2, 0, 2, 0)
        painter.drawLine(0, -2, 0, 2)
    
                 
class MapperScene(QtWidgets.QGraphicsScene):

    def __init__(self, **kwds):
        super().__init__(**kwds)
        self.channels = []

    def updateDisplay(self, current_channel, channels, groups, mappings, frame_no, fmin, fmax):
        self.clear()

        current_channel.addFrame(self, frame_no, fmin, fmax)

        for channel in channels:
            channel.addLocs(self, mappings)

        for group in groups:
            group.addGroup(self)
                 
    
class Window(QtWidgets.QMainWindow):

    def __init__(self, **kwds):
        super().__init__(**kwds)

        self.colors = [QtGui.QColor(255,0,0),
                       QtGui.QColor(0,255,0),
                       QtGui.QColor(0,0,255),
                       QtGui.QColor(255,255,0),
                       QtGui.QColor(255,0,255),
                       QtGui.QColor(0,255,255)]
        self.channels = []
        self.current_channel = None
        self.directory = ""
        self.groups = []
        self.mappings = {}
        self.settings = QtCore.QSettings("Zhuang Lab", "mapper")

        self.ui = mapperUi.Ui_MainWindow()
        self.ui.setupUi(self)

        # Setup graphics view.
        self.mapper_scene = MapperScene()
        self.ui.channelGraphicsView.setScene(self.mapper_scene)

        # Load settings.
        self.directory = str(self.settings.value("directory", ""))
        self.move(self.settings.value("position", self.pos()))
        self.resize(self.settings.value("size", self.size()))
        self.ui.maxDistDoubleSpinBox.setValue(float(self.settings.value("max_dist", 8.0)))
        self.ui.maxSpinBox.setValue(int(self.settings.value("maximum", 2000)))
        self.ui.minSpinBox.setValue(int(self.settings.value("minimum", 100)))
        
        # Add range slider.
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
        
        # Connect signals.
        self.ui.actionClear_Groups.triggered.connect(self.handleClearGroups)
        self.ui.actionLoad_Channel.triggered.connect(self.handleLoadChannel)
        self.ui.actionReset.triggered.connect(self.handleReset)
        self.ui.actionSave_Mapping.triggered.connect(self.handleSaveMapping)
        self.ui.actionQuit.triggered.connect(self.handleQuit)
        self.ui.channelComboBox.currentIndexChanged.connect(self.handleChannelComboBox)
        self.ui.channelGraphicsView.mouseClick.connect(self.handleChannelGraphicsView)
        self.ui.flipLRPushButton.clicked.connect(self.handleFlipLR)
        self.ui.flipUDPushButton.clicked.connect(self.handleFlipUD)
        self.ui.mapItPushButton.clicked.connect(self.handleMapIt)
        self.ui.maxSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.minSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)

    def cleanUp(self):
        self.settings.setValue("directory", self.directory)
        self.settings.setValue("max_dist", self.ui.maxDistDoubleSpinBox.value())
        self.settings.setValue("maximum", self.ui.maxSpinBox.value())
        self.settings.setValue("minimum", self.ui.minSpinBox.value())
        self.settings.setValue("position", self.pos())
        self.settings.setValue("size", self.size())

    def closeEvent(self, event):
        self.cleanUp()

    def handleChannelComboBox(self, index):
        if self.current_channel is None:
            return

        self.current_channel = self.channels[index]
        self.updateScene()

    def handleChannelGraphicsView(self, click_x, click_y):
        if (len(self.channels) < 2):
            return
        
        group = Group()
        for channel in self.channels:
            coords = channel.findNearest(click_x,
                                         click_y,
                                         self.ui.maxDistDoubleSpinBox.value())
            if coords is not None:
                group.addLocalization(coords[0], coords[1])

        # Check that the group has a representative from each channel.
        if (group.getSize() == len(self.channels)):

            # Check that we don't already have this group.
            for a_group in self.groups:
                if group.isSame(a_group):
                    return
            
            self.groups.append(group)

            self.updateScene()

    def handleClearGroups(self):
        self.groups = []
        self.mappings = {}
        self.updateScene()

    def handleFlipLR(self):
        if self.current_channel is not None:
            if (self.current_channel != self.channels[0]):
                self.current_channel.flipLR()
                self.updateScene()

    def handleFlipUD(self):
        if self.current_channel is not None:
            if (self.current_channel != self.channels[0]):
                self.current_channel.flipUD()
                self.updateScene()

    def handleLoadChannel(self):
        movie_name = QtWidgets.QFileDialog.getOpenFileName(self,
                                                           "Channel Movie",
                                                           self.directory,
                                                           "*.dax *.spe *.tif")[0]        
        if (len(movie_name) == 0):
            return
        
        self.directory = os.path.dirname(movie_name)
        
        locs_name = QtWidgets.QFileDialog.getOpenFileName(self,
                                                          "Load Localizations",
                                                          self.directory,
                                                          "*.bin *.hdf5")[0]
        if (len(locs_name) == 0):
            return

        self.directory = os.path.dirname(locs_name)

        index = len(self.channels)
        channel = Channel(color = self.colors[index],
                          locs_name = locs_name,
                          movie_name = movie_name,
                          number = index)
        self.current_channel = channel
        self.channels.append(channel)

        self.groups = []

        self.ui.channelComboBox.addItem("Channel" + str(index))
        self.ui.channelComboBox.setCurrentIndex(index)

    def handleMapIt(self):
        if (len(self.groups) < 2):
            return

        self.mappings = {}

        # 0 <-> 0 is the identity matrix.
        self.mappings["0_0_x"] = numpy.array([0.0,1.0,0.0])
        self.mappings["0_0_y"] = numpy.array([0.0,0.0,1.0])

        # Helper function creating transform matrix and vectors.
        def getXYM(ch1, ch2):
            x = numpy.zeros(len(self.groups))
            y = numpy.zeros(len(self.groups))
            m = numpy.ones((len(self.groups), 3))
            for j, group in enumerate(self.groups):
                [ch1_x, ch1_y] = group.getXY(ch1)
                [ch2_x, ch2_y] = group.getXY(ch2)
                x[j] = ch1_x
                y[j] = ch1_y
                m[j,1] = ch2_x
                m[j,2] = ch2_y
            return [x, y, m]

        # Calculate transform matrix for each channel to channel 0.
        for i in range(len(self.channels) - 1):

            # Channel i -> Channel 0
            [x, y, m] = getXYM(0, i + 1)
            m_name = str(i+1) + "_" + str(0)
            self.mappings[m_name + "_x"] = numpy.linalg.lstsq(m, x)[0]
            self.mappings[m_name + "_y"] = numpy.linalg.lstsq(m, y)[0]
            
            # Channel 0 -> Channel i
            [x, y, m] = getXYM(i + 1, 0)
            m_name = str(0) + "_" + str(i+1)
            self.mappings[m_name + "_x"] = numpy.linalg.lstsq(m, x)[0]
            self.mappings[m_name + "_y"] = numpy.linalg.lstsq(m, y)[0]

        #
        # The mapping name convention is 'from channel' underscore 'to channel',
        # i.e. the mapping 1_0_x will transform the x coordinate of channel1 into
        # channel0 coordinates. Mapping are all of the form:
        #
        # xf = cx1 + cx2 * xi + cx3 * yi.
        # yf = cy1 + cy2 * xi + cx3 * yi.
        #
        for key in sorted(self.mappings):
            print(key, self.mappings[key])
            
        self.updateScene()
            
    def handleMaxMinSpinBox(self, value):
        self.rangeSlider.setValues([self.ui.minSpinBox.value(),
                                    self.ui.maxSpinBox.value()])
        
    def handleQuit(self):
        self.close()

    def handleRangeChange(self, range_min, range_max):
        self.ui.minSpinBox.setValue(range_min)
        self.ui.maxSpinBox.setValue(range_max)
        self.updateScene()

    def handleReset(self):
        self.current_channel = None
        self.channels = []
        self.groups = []
        self.mapper_scene.clear()
        self.mappings = {}
        self.ui.frameLabel.setText("")

    def handleSaveMapping(self):
        if (len(self.mappings) == 0):
            return

        mapping_name = QtWidgets.QFileDialog.getSaveFileName(self,
                                                             "Channel Mappings",
                                                             self.directory,
                                                             "*.map")[0]
        if (len(mapping_name) > 0):
            with open(mapping_name, 'wb') as fp:
                pickle.dump(self.mappings, fp)

    def keyPressEvent(self, event):
        if self.current_channel is None:
            return
        
        cf = self.current_channel.getCurrentFrameNumber()

        # These go back and forth through the frames.
        if (event.key() == QtCore.Qt.Key_End):
            self.current_channel.changeFrame(self.current_channel.getMovieLength())
            self.groups = []
        elif (event.key() == QtCore.Qt.Key_Home):
            self.current_channel.changeFrame(0)
            self.groups = []
        elif (event.key() == QtCore.Qt.Key_Comma):
            self.current_channel.changeFrame(cf - 1)
            self.groups = []
        elif (event.key() == QtCore.Qt.Key_Period):
            self.current_channel.changeFrame(cf + 1)
            self.groups = []
        elif (event.key() == QtCore.Qt.Key_K):
            self.current_channel.changeFrame(cf - 20)
            self.groups = []
        elif (event.key() == QtCore.Qt.Key_L):
            self.current_channel.changeFrame(cf + 20)
            self.groups = []

        # These change the channel X/Y offset. Channel 0 cannot be offset.
        if (self.current_channel != self.channels[0]):
            if (event.key() == QtCore.Qt.Key_A):
                self.current_channel.offsetX(-1)
                self.groups = []
            elif (event.key() == QtCore.Qt.Key_D):
                self.current_channel.offsetX(1)
                self.groups = []
            elif (event.key() == QtCore.Qt.Key_W):
                self.current_channel.offsetY(-1)
                self.groups = []
            elif (event.key() == QtCore.Qt.Key_S):
                self.current_channel.offsetY(1)
                self.groups = []

        self.updateScene()

    def updateScene(self):
        cf = self.current_channel.getCurrentFrameNumber()
        ml = self.current_channel.getMovieLength()
        self.ui.frameLabel.setText("Frame: " + str(cf + 1) + " / " + str(ml))
        self.ui.groupsLabel.setText(str(len(self.groups)) + " groups")
        self.mapper_scene.updateDisplay(self.current_channel,
                                        self.channels,
                                        self.groups,
                                        self.mappings,
                                        cf,
                                        self.ui.minSpinBox.value(),
                                        self.ui.maxSpinBox.value())


if (__name__ == "__main__"):
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()
    
#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
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
