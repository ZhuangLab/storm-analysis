#!/usr/bin/env python
"""
QGraphics specialized for mapper.

Hazen 05/17
"""

from PyQt5 import QtCore, QtGui, QtWidgets


class MapperView(QtWidgets.QGraphicsView):
    mouseClick = QtCore.pyqtSignal(float, float)

    def __init__(self, parent):
        super().__init__(parent)
        self.scale_int = 0

        self.setRenderHint(QtGui.QPainter.Antialiasing + QtGui.QPainter.SmoothPixmapTransform)
        self.scale(1.0, 1.0)
        
    def mousePressEvent(self, event):
        pointf = self.mapToScene(event.pos())
        self.mouseClick.emit(pointf.x(), pointf.y())

    def wheelEvent(self, event):
        if not event.angleDelta().isNull():
            if (event.angleDelta().y() > 0):
                self.scale_int += 1
            else:
                self.scale_int -= 1

            if (self.scale_int == 0):
                flt_scale = 1.0
            elif (self.scale_int > 0):
                flt_scale = self.scale_int + 1.0
            else:
                flt_scale = 1.0/(-self.scale_int + 1.0)

            transform = QtGui.QTransform()
            transform.scale(flt_scale, flt_scale)
            self.setTransform(transform)
            event.accept()
