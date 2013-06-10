#!/usr/bin/python
#
# Qt Widget Range slider widget.
#
# Hazen 06/13
#

from PyQt4 import QtCore, QtGui
import sys

# Range Slider super class
class QRangeSlider(QtGui.QWidget):
    doubleClick = QtCore.pyqtSignal()
    rangeChanged = QtCore.pyqtSignal(float, float)

    def __init__(self, slider_range, values, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.bar_width = 10
        self.emit_while_moving = 0
        self.moving = "none"
        self.old_scale_min = 0.0
        self.old_scale_max = 0.0
        self.scale = 0
        self.setMouseTracking(False)
        self.single_step = 0.0

        if slider_range:
            self.setRange(slider_range)
        else:
            self.setRange([0.0, 1.0, 0.01])
        if values:
            self.setValues(values)
        else:
            self.setValues([0.3, 0.6])

        self.setFocusPolicy(QtCore.Qt.ClickFocus)

    def emitRange(self):
        if (self.old_scale_min != self.scale_min) or (self.old_scale_max != self.scale_max):
            self.rangeChanged.emit(self.scale_min, self.scale_max)
            self.old_scale_min = self.scale_min
            self.old_scale_max = self.scale_max
            if 0:
                print "Range change:", self.scale_min, self.scale_max

    def getValues(self):
        return [self.scale_min, self.scale_max]

    def keyPressEvent(self, event):
        key = event.key()

        # move bars based on arrow keys
        moving_max = False
        if key == QtCore.Qt.Key_Up:
            self.scale_max += self.single_step
            moving_max = True
        elif key == QtCore.Qt.Key_Down:
            self.scale_max -= self.single_step
            moving_max = True
        elif key == QtCore.Qt.Key_Left:
            self.scale_min -= self.single_step
        elif key == QtCore.Qt.Key_Right:
            self.scale_min += self.single_step

        # update (if necessary) based on allowed range
        if moving_max:
            if (self.scale_max < self.scale_min):
                self.scale_min = self.scale_max
        else:
            if (self.scale_min > self.scale_max):
                self.scale_max = self.scale_min

        if (self.scale_min < self.start):
            self.scale_min = self.start
        if (self.scale_max < self.start):
            self.scale_max = self.start

        slider_max = self.start + self.scale
        if (self.scale_min > slider_max):
            self.scale_min = slider_max
        if (self.scale_max > slider_max):
            self.scale_max = slider_max

        self.emitRange()
        self.updateDisplayValues()
        self.update()
    
    def mouseDoubleClickEvent(self, event):
        #self.emit(QtCore.SIGNAL("doubleClick()"))
        self.doubleClick.emit()

    def mouseMoveEvent(self, event):
        size = self.rangeSliderSize()
        diff = self.start_pos - self.getPos(event)
        if self.moving == "min":
            temp = self.start_display_min - diff
            if (temp >= self.bar_width) and (temp < size - self.bar_width):
                self.display_min = temp
                if self.display_max < self.display_min:
                    self.display_max = self.display_min
                self.updateScaleValues()
                if self.emit_while_moving:
                    self.emitRange()
        elif self.moving == "max":
            temp = self.start_display_max - diff
            if (temp >= self.bar_width) and (temp < size - self.bar_width):
                self.display_max = temp
                if self.display_max < self.display_min:
                    self.display_min = self.display_max
                self.updateScaleValues()
                if self.emit_while_moving:
                    self.emitRange()
        elif self.moving == "bar":
            temp = self.start_display_min - diff
            if (temp >= self.bar_width) and (temp < size - self.bar_width - (self.start_display_max - self.start_display_min)):
                self.display_min = temp
                self.display_max = self.start_display_max - diff
                self.updateScaleValues()
                if self.emit_while_moving:
                    self.emitRange()

    def mousePressEvent(self, event):
        pos = self.getPos(event)
        if abs(self.display_min - 0.5 * self.bar_width - pos) < (0.5 * self.bar_width):
            self.moving = "min"
        elif abs(self.display_max + 0.5 * self.bar_width - pos) < (0.5 * self.bar_width):
            self.moving = "max"
        elif (pos > self.display_min) and (pos < self.display_max):
            self.moving = "bar"
        self.start_display_min = self.display_min
        self.start_display_max = self.display_max
        self.start_pos = pos
            
    def mouseReleaseEvent(self, event):
        if not (self.moving == "none"):
            self.emitRange()
        self.moving = "none"

    def resizeEvent(self, event):
        self.updateDisplayValues()

    def setRange(self, slider_range):
        self.start = slider_range[0]
        self.scale = slider_range[1] - slider_range[0]
        self.single_step = slider_range[2]

    def setValues(self, values):
        self.scale_min = values[0]
        self.scale_max = values[1]
        self.emitRange()
        self.updateDisplayValues()
        self.update()

    def setEmitWhileMoving(self, flag):
        if flag:
            self.emit_while_moving = 1
        else:
            self.emit_while_moving = 0

    def updateDisplayValues(self):
        size = float(self.rangeSliderSize() - 2 * self.bar_width - 1)
        self.display_min = int(size * (self.scale_min - self.start)/self.scale) + self.bar_width
        self.display_max = int(size * (self.scale_max - self.start)/self.scale) + self.bar_width

    def updateScaleValues(self):
        size = float(self.rangeSliderSize() - 2 * self.bar_width - 1)
        if (self.moving == "min") or (self.moving == "bar"):
            self.scale_min = self.start + (self.display_min - self.bar_width)/float(size) * self.scale
            self.scale_min = float(round(self.scale_min/self.single_step))*self.single_step
        if (self.moving == "max") or (self.moving == "bar"):
            self.scale_max = self.start + (self.display_max - self.bar_width)/float(size) * self.scale
            self.scale_max = float(round(self.scale_max/self.single_step))*self.single_step
        self.updateDisplayValues()
        self.update()

# Horizontal Range Slider
class QHRangeSlider(QRangeSlider):
    def __init__(self, slider_range = None, values = None, parent = None):
        QRangeSlider.__init__(self, slider_range, values, parent)
        if (not parent):
            self.setGeometry(200, 200, 200, 20)

    def getPos(self, event):
        return event.x()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        w = self.width()
        h = self.height()

        # background
        painter.setPen(QtCore.Qt.gray)
        painter.setBrush(QtCore.Qt.lightGray)
        painter.drawRect(2, 2, w-4, h-4)

        # range bar
        painter.setPen(QtCore.Qt.darkGray)
        painter.setBrush(QtCore.Qt.darkGray)
        painter.drawRect(self.display_min-1, 5, self.display_max-self.display_min+2, h-10)

        # min & max tabs
        painter.setPen(QtCore.Qt.black)
        painter.setBrush(QtCore.Qt.gray)
        painter.drawRect(self.display_min-self.bar_width, 1, self.bar_width, h-2)

        painter.setPen(QtCore.Qt.black)
        painter.setBrush(QtCore.Qt.gray)
        painter.drawRect(self.display_max, 1, self.bar_width, h-2)

    def rangeSliderSize(self):
        return self.width()

# Vertical Range Slider
class QVRangeSlider(QRangeSlider):
    def __init__(self, slider_range = None, values = None, parent = None):
        QRangeSlider.__init__(self, slider_range, values, parent)
        if (not parent):
            self.setGeometry(200, 200, 20, 200)

    def getPos(self, event):
        return self.height() - event.y()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        w = self.width()
        h = self.height()

        # background
        painter.setPen(QtCore.Qt.gray)
        painter.setBrush(QtCore.Qt.lightGray)
        painter.drawRect(2, 2, w-4, h-4)

        # range bar
        painter.setPen(QtCore.Qt.darkGray)
        painter.setBrush(QtCore.Qt.darkGray)
        painter.drawRect(5, h-self.display_max-1, w-10, self.display_max-self.display_min+1)

        # min & max tabs
        painter.setPen(QtCore.Qt.black)
        painter.setBrush(QtCore.Qt.gray)
        painter.drawRect(1, h-self.display_max-self.bar_width-1, w-2, self.bar_width)

        painter.setPen(QtCore.Qt.black)
        painter.setBrush(QtCore.Qt.gray)
        painter.drawRect(1, h-self.display_min-1, w-2, self.bar_width)

    def rangeSliderSize(self):
        return self.height()


#
# Testing
#

if __name__ == "__main__":
    class Parameters:
        def __init__(self):
            self.x_pixels = 200
            self.y_pixels = 200

    app = QtGui.QApplication(sys.argv)
    if 0:
        hslider = QHRangeSlider(slider_range = [-5.0, 5.0, 0.5], values = [-2.5, 2.5])
        hslider.setEmitWhileMoving(True)
        hslider.show()
    if 1:
        vslider = QVRangeSlider(slider_range = [-5.0, 5.0, 0.5], values = [-2.5, 2.5])
        vslider.setEmitWhileMoving(True)
        vslider.show()
    sys.exit(app.exec_())


#
# The MIT License
#
# Copyright (c) 2009 Zhuang Lab, Harvard University
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

