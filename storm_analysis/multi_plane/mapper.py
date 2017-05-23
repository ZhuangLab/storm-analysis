#!/usr/bin/env python
"""
Utility for identifying the transform between different
channels in a multi-channel data set.

This assumes that each channel is it's own movie and
that the localizations in the channel have been 
identified with a tool like 3D-DAOSTORM.

Hazen 05/17
"""

import numpy
import os
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.visualizer.qtRangeSlider as qtRangeSlider
import storm_analysis.sa_library.readinsight3 as readinsight3


import storm_analysis.multi_plane.mapper_ui as mapperUi


class Channel(object):
    """
    Wrapper for a STORM movie and it's associated localization list.
    """
    def __init__(self, color = None, locs_name = None, movie_name = None, **kwds):
        super().__init__(**kwds)
        self.color = color
        self.cur_frame = 0
        self.flip_lr = False
        self.flip_ud = False
        self.fr_height = None
        self.fr_width = None
        self.locs_name = locs_name
        self.movie_name = movie_name
        self.offset_x = 0
        self.offset_y = 0

        self.locs_i3 = readinsight3.I3Reader(locs_name)
        self.movie_fp = datareader.inferReader(movie_name)
        self.movie_len = self.movie_fp.filmSize()[2]

    def addFrame(self, scene, frame_no, fmin, fmax):
        """
        Add a frame to a graphics scene.
        """
        frame = self.movie_fp.loadAFrame(self.cur_frame).astype(numpy.float)
        frame = numpy.transpose(frame)

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

        image = QtGui.QImage(frame_RGB.data,
                             self.fr_width,
                             self.fr_height,
                             QtGui.QImage.Format_RGB32)
        image.ndarray1 = frame
        image.ndarray2 = frame_RGB
    
        # Add to scene
        pixmap = QtGui.QPixmap.fromImage(image.mirrored(self.flip_lr, self.flip_ud))
        pixmap_item = QtWidgets.QGraphicsPixmapItem(pixmap)
        pixmap_item.setOffset(self.offset_x, self.offset_y)
        scene.addItem(pixmap_item)

    def addLocs(self, scene):
        [x, y] = self.getLocs()
        
        for i in range(x.size):
            loc_item = LocalizationItem(i, x[i], y[i], 6.0, self.color)
            scene.addItem(loc_item)

    def changeFrame(self, frame_no):
        if (frame_no < 0):
            self.cur_frame = 0
        elif (frame_no >= self.movie_len):
            self.cur_frame = self.movie_len - 1
        else:
            self.cur_frame = frame_no
        return self.cur_frame

    def flipLR(self):
        self.flip_lr = not self.flip_lr

    def flipUD(self):
        self.flip_ud = not self.flip_ud
        
    def getCurrentFrameNumber(self):
        return self.cur_frame

    def getLocs(self):
        locs = self.locs_i3.getMoleculesInFrame(self.cur_frame+1)
        x = locs["x"]
        y = locs["y"]

        if self.flip_lr:
            x = self.fr_width - x + 1 + self.offset_x
        else:
            x += self.offset_x
            
        if self.flip_ud:
            y = self.fr_height - y + 1 + self.offset_y
        else:
            y += self.offset_y

        return [x, y]
        
    def getMovieLength(self):
        return self.movie_len

    def offsetX(self, inc):
        self.offset_x += inc

    def offsetY(self, inc):
        self.offset_y += inc


class LocalizationItem(QtWidgets.QGraphicsEllipseItem):
    """
    Localization item for a graphics scene.
    """
    def __init__(self, loc_index, x, y, d, color, **kwds):
        self.loc_index = loc_index

        self.default_pen = QtGui.QPen(color)        
        self.default_pen.setWidthF(0.6)

        x = x - 0.5*d - 0.5
        y = y - 0.5*d - 0.5
        QtWidgets.QGraphicsEllipseItem.__init__(self, x, y, d, d)
        self.setPen(self.default_pen)


class MapperScene(QtWidgets.QGraphicsScene):

    def __init__(self, **kwds):
        super().__init__(**kwds)
        self.channels = []

    def updateDisplay(self, current_channel, channels, frame_no, fmin, fmax):
        self.clear()

        current_channel.addFrame(self, frame_no, fmin, fmax)

        for channel in channels:
            channel.addLocs(self)
                 
    
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
        self.ui.maxSpinBox.setValue(int(self.settings.value("maximum", 2000)))
        self.ui.minSpinBox.setValue(int(self.settings.value("minimum", 100)))
        self.ui.toleranceDoubleSpinBox.setValue(float(self.settings.value("tolerance", 1.0)))
        
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
        self.ui.actionLoad_Channel.triggered.connect(self.handleLoadChannel)
        self.ui.actionReset.triggered.connect(self.handleReset)
        self.ui.actionQuit.triggered.connect(self.handleQuit)
        self.ui.channelComboBox.currentIndexChanged.connect(self.handleChannelComboBox)
        self.ui.channelGraphicsView.mouseClick.connect(self.handleChannelGraphicsView)
        self.ui.flipLRPushButton.clicked.connect(self.handleFlipLR)
        self.ui.flipUDPushButton.clicked.connect(self.handleFlipUD)
        self.ui.maxSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)
        self.ui.minSpinBox.valueChanged.connect(self.handleMaxMinSpinBox)

    def cleanUp(self):
        self.settings.setValue("directory", self.directory)
        self.settings.setValue("maximum", self.ui.maxSpinBox.value())
        self.settings.setValue("minimum", self.ui.minSpinBox.value())
        self.settings.setValue("position", self.pos())
        self.settings.setValue("size", self.size())
        self.settings.setValue("tolerance", self.ui.toleranceDoubleSpinBox.value())

    def closeEvent(self, event):
        self.cleanUp()

    def handleChannelComboBox(self, index):
        if self.current_channel is None:
            return

        self.current_channel = self.channels[index]
        self.updateScene()

    def handleChannelGraphicsView(self, click_x, click_y):
        print(click_x, click_y)

    def handleFlipLR(self):
        if self.current_channel is not None:
            self.current_channel.flipLR()
            self.updateScene()

    def handleFlipUD(self):
        if self.current_channel is not None:
            self.current_channel.flipUD()
            self.updateScene()

    def handleLoadChannel(self):
        movie_name = QtWidgets.QFileDialog.getOpenFileName(self,
                                                           "Channel Movie",
                                                           self.directory,
                                                           "*.dax *.spe *.tif")[0]
        if movie_name is None:
            return
        self.directory = os.path.dirname(movie_name)
        
        locs_name = QtWidgets.QFileDialog.getOpenFileName(self,
                                                          "Load Localizations",
                                                          self.directory,
                                                          "*.bin")[0]        
        if locs_name is None:
            return
        self.directory = os.path.dirname(locs_name)

        index = len(self.channels)
        channel = Channel(color = self.colors[index],
                          locs_name = locs_name,
                          movie_name = movie_name)
        self.current_channel = channel
        self.channels.append(channel)

        self.ui.channelComboBox.addItem("Channel" + str(index))
        self.ui.channelComboBox.setCurrentIndex(index)
#        self.updateScene()
        
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
        self.mapper_scene.clear()
        self.ui.frameLabel.setText("")

    def keyPressEvent(self, event):
        if self.current_channel is None:
            return
        
        cf = self.current_channel.getCurrentFrameNumber()

        # These go back and forth through the frames.
        if (event.key() == QtCore.Qt.Key_End):
            self.current_channel.changeFrame(self.current_channel.getMovieLength())
        elif (event.key() == QtCore.Qt.Key_Home):
            self.current_channel.changeFrame(0)
        elif (event.key() == QtCore.Qt.Key_Comma):
            self.current_channel.changeFrame(cf - 1)
        elif (event.key() == QtCore.Qt.Key_Period):
            self.current_channel.changeFrame(cf + 1)
        elif (event.key() == QtCore.Qt.Key_K):
            self.current_channel.changeFrame(cf - 20)
        elif (event.key() == QtCore.Qt.Key_L):
            self.current_channel.changeFrame(cf + 20)

        # These change the channel X/Y offset.
        elif (event.key() == QtCore.Qt.Key_A):
            self.current_channel.offsetX(-1)
        elif (event.key() == QtCore.Qt.Key_D):
            self.current_channel.offsetX(1)
        elif (event.key() == QtCore.Qt.Key_W):
            self.current_channel.offsetY(-1)
        elif (event.key() == QtCore.Qt.Key_S):
            self.current_channel.offsetY(1)

        self.updateScene()

    def updateScene(self):
        cf = self.current_channel.getCurrentFrameNumber()
        ml = self.current_channel.getMovieLength()
        self.ui.frameLabel.setText("Frame: " + str(cf + 1) + " / " + str(ml))
        self.mapper_scene.updateDisplay(self.current_channel,
                                        self.channels,
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
