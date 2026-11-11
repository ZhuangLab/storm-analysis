#!/usr/bin/env python
"""
Tests of the visualizer GUI. These need PyQt6 and are skipped without it.

They run on Qt's offscreen platform, so nothing is drawn on the screen.
"""
import os
import sys
import warnings

import numpy
import pytest

pytest.importorskip("PyQt6.QtWidgets")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6 import QtCore, QtGui, QtTest, QtWidgets

import storm_analysis
import storm_analysis.sa_library.datawriter as datawriter
import storm_analysis.sa_library.sa_h5py as saH5Py


# Qt needs the application to outlive every widget, so keep it here rather
# than in a local variable that is garbage collected.
application = None

def getApplication():
    global application
    if application is None:
        application = QtWidgets.QApplication.instance()
    if application is None:
        application = QtWidgets.QApplication([])
    return application


class RecordExceptions(object):
    """
    PyQt calls sys.excepthook for an exception raised inside a slot or an
    event handler, and aborts the process if the hook is the default one.
    This records them instead so that a test can fail on them.
    """
    def __enter__(self):
        self.exceptions = []
        self.original_hook = sys.excepthook
        sys.excepthook = lambda exception_type, exception, traceback: self.exceptions.append(exception)
        return self.exceptions

    def __exit__(self, exception_type, exception, traceback):
        sys.excepthook = self.original_hook


def makeWindow():
    """
    Create a visualizer window that starts from default settings, and
    reads and writes them in the test output directory rather than in
    the user's own settings.

    Note that QSettings.setDefaultFormat() and setPath() can't do this,
    as QSettings("Zhuang Lab", "visualizer") always uses the native
    format.
    """
    getApplication()
    settings_name = storm_analysis.getPathOutputTest("test_visualizer_settings.ini")
    storm_analysis.removeFile(settings_name)
    settings = QtCore.QSettings(settings_name, QtCore.QSettings.Format.IniFormat)

    import storm_analysis.visualizer.visualizer as visualizer
    window = visualizer.Window(settings = settings)
    window.reported_errors = []
    window.showError = window.reported_errors.append
    window.show()
    QtWidgets.QApplication.processEvents()
    return window


def makeTestData(name, localizations, frame_shape = (32, 32), n_frames = 2):
    """
    Write a movie and a matching HDF5 localization file, with the
    localizations in the first frame. Returns the two file names.
    """
    movie_name = storm_analysis.getPathOutputTest(name + ".dax")
    h5_name = storm_analysis.getPathOutputTest(name + ".hdf5")
    storm_analysis.removeFile(h5_name)

    frame = numpy.full(frame_shape, 100.0)
    frame[2, 3] = 900.0
    movie = datawriter.DaxWriter(movie_name)
    for i in range(n_frames):
        movie.addFrame(frame)
    movie.close()

    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(frame_shape[1], frame_shape[0], n_frames, "")
        h5.addLocalizations(localizations, 0)

    return movie_name, h5_name


def showLocalizationsNow(window):
    """
    Show the localizations without waiting for the timer that delays them
    after a frame change.
    """
    window.locs_display_timer.stop()
    window.handleLocsDisplayTimer()


def tableValue(table_widget, field):
    for row in range(table_widget.rowCount()):
        if (table_widget.item(row, 0).text() == field):
            return table_widget.item(row, 1).text()
    raise KeyError(field)


def clickAt(movie_view, scene_x, scene_y):
    position = movie_view.mapFromScene(QtCore.QPointF(scene_x, scene_y))
    QtTest.QTest.mouseClick(movie_view.viewport(),
                            QtCore.Qt.MouseButton.LeftButton,
                            QtCore.Qt.KeyboardModifier.NoModifier,
                            position)


def twoNeighbors():
    """
    Two localizations one pixel apart on the diagonal.
    """
    return {"x" : numpy.array([10.0, 11.0]),
            "y" : numpy.array([10.0, 11.0]),
            "xsigma" : numpy.array([1.0, 1.0]),
            "height" : numpy.array([100.0, 200.0])}


def test_visualizer_load_and_step():
    """
    Load a movie and its localizations, then step through the frames.
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    h5_name = storm_analysis.getData("test/data/test_ref.hdf5")

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.loadLocalizations(h5_name, 1)
            showLocalizationsNow(window)
            assert (window.ui.frameLabel.text() == "frame 1 (10)")

            with saH5Py.SAH5Py(h5_name) as h5:
                n_locs = h5.getLocalizationsInFrame(0)["x"].size
            assert (n_locs > 0)
            assert (len(window.locs1_list.mol_items) == n_locs)
            for item in window.locs1_list.mol_items:
                assert (item.scene() is window.movie_view.scene)

            QtTest.QTest.keyClick(window.movie_view, QtCore.Qt.Key.Key_Period)
            assert (window.ui.frameLabel.text() == "frame 2 (10)")
            QtTest.QTest.keyClick(window.movie_view, QtCore.Qt.Key.Key_End)
            assert (window.ui.frameLabel.text() == "frame 10 (10)")
            QtTest.QTest.keyClick(window.movie_view, QtCore.Qt.Key.Key_Home)
            assert (window.ui.frameLabel.text() == "frame 1 (10)")
        finally:
            window.close()

    assert (exceptions == [])
    assert (window.reported_errors == [])


def test_visualizer_click_selects_localization():
    """
    Clicking on a localization's circle shows that localization's fit.
    """
    movie_name, h5_name = makeTestData("test_visualizer_click", twoNeighbors())

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.loadLocalizations(h5_name, 1)
            showLocalizationsNow(window)
            movie_view = window.movie_view
            movie_view.setTransform(QtGui.QTransform.fromScale(10.0, 10.0))
            movie_view.centerOn(16.0, 16.0)

            # The circle of a localization at (x, y) is centered on pixel (x, y).
            for item, x in zip(window.locs1_list.mol_items, [10.0, 11.0]):
                center = item.rect().center()
                assert (abs(center.x() - (x + 0.5)) < 1.0e-6)
                assert (abs(center.y() - (x + 0.5)) < 1.0e-6)

            clickAt(movie_view, 10.5, 10.5)
            assert (tableValue(window.ui.locs1TableWidget, "x") == "10.000")
            assert (tableValue(window.ui.locs1TableWidget, "height") == "100.000")

            clickAt(movie_view, 11.5, 11.5)
            assert (tableValue(window.ui.locs1TableWidget, "x") == "11.000")
            assert (tableValue(window.ui.locs1TableWidget, "height") == "200.000")

            # The position readout uses the same convention.
            QtTest.QTest.mouseMove(movie_view.viewport(),
                                   movie_view.mapFromScene(QtCore.QPointF(3.5, 2.5)))
            assert (window.ui.xyiLabel.text() == "3.00, 2.00, 900")
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_click_after_contrast_change():
    """
    Changing the contrast keeps the localizations, so a click right after
    a contrast change still finds one.
    """
    movie_name, h5_name = makeTestData("test_visualizer_contrast", twoNeighbors())

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.loadLocalizations(h5_name, 1)
            showLocalizationsNow(window)
            window.movie_view.setTransform(QtGui.QTransform.fromScale(10.0, 10.0))
            window.movie_view.centerOn(16.0, 16.0)

            window.ui.maxSpinBox.setValue(1500)
            for item in window.locs1_list.mol_items:
                assert (item.scene() is window.movie_view.scene)

            clickAt(window.movie_view, 10.5, 10.5)
            assert (tableValue(window.ui.locs1TableWidget, "x") == "10.000")
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_equal_contrast_limits():
    """
    Setting the minimum equal to the maximum thresholds the image instead
    of dividing by zero.
    """
    movie_name, h5_name = makeTestData("test_visualizer_threshold", twoNeighbors())

    with RecordExceptions() as exceptions, warnings.catch_warnings(record = True) as caught_warnings:
        warnings.simplefilter("always")
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.ui.minSpinBox.setValue(500)
            window.ui.maxSpinBox.setValue(500)
            QtWidgets.QApplication.processEvents()
            image = window.movie_view.pixmap_item.pixmap().toImage()
            assert (image.pixelColor(3, 2).value() == 255)
            assert (image.pixelColor(0, 0).value() == 0)
        finally:
            window.close()

    assert (exceptions == [])
    assert ([str(warning.message) for warning in caught_warnings if issubclass(warning.category, RuntimeWarning)] == [])


def test_visualizer_bad_files_are_reported():
    """
    A file that can't be read is reported and doesn't stop the program.
    """
    missing_movie = storm_analysis.getPathOutputTest("test_visualizer_missing.dax")
    missing_locs = storm_analysis.getPathOutputTest("test_visualizer_missing.hdf5")
    unsupported = storm_analysis.getPathOutputTest("test_visualizer_unsupported.txt")
    for name in [missing_movie, missing_locs]:
        storm_analysis.removeFile(name)
    with open(unsupported, "w") as text_file:
        text_file.write("not a movie\n")

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(unsupported)
            window.loadMovie(missing_movie)
            assert (window.movie_file is None)
            assert (len(window.reported_errors) == 2)

            movie_name, h5_name = makeTestData("test_visualizer_bad_files", twoNeighbors())
            window.loadMovie(movie_name)
            window.loadLocalizations(missing_locs, 2)
            assert (window.locs2_list is None)
            assert (len(window.reported_errors) == 3)
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_new_movie_closes_old_one():
    """
    Loading a movie closes the file of the movie it replaces.
    """
    movie_name_1, h5_name_1 = makeTestData("test_visualizer_first", twoNeighbors())
    movie_name_2, h5_name_2 = makeTestData("test_visualizer_second", twoNeighbors())

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name_1)
            first_movie = window.movie_file
            window.loadMovie(movie_name_2)
            assert (first_movie.fileptr is None)
            second_movie = window.movie_file
        finally:
            window.close()
        assert (second_movie.fileptr is None)

    assert (exceptions == [])


def test_visualizer_integer_fields():
    """
    Integer fields are shown as integers whatever their width.
    """
    localizations = twoNeighbors()
    localizations["category"] = numpy.array([3, 4], dtype = numpy.int32)
    localizations["iterations"] = numpy.array([7, 8], dtype = numpy.int64)
    movie_name, h5_name = makeTestData("test_visualizer_integers", localizations)

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.loadLocalizations(h5_name, 1)
            showLocalizationsNow(window)
            window.movie_view.setTransform(QtGui.QTransform.fromScale(10.0, 10.0))
            window.movie_view.centerOn(16.0, 16.0)

            clickAt(window.movie_view, 10.5, 10.5)
            assert (tableValue(window.ui.locs1TableWidget, "category") == "3")
            assert (tableValue(window.ui.locs1TableWidget, "iterations") == "7")
            assert (tableValue(window.ui.locs1TableWidget, "x") == "10.000")
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_movie_fits_view():
    """
    A newly loaded movie is scaled to fit the view.
    """
    movie_name, h5_name = makeTestData("test_visualizer_fit", twoNeighbors(), frame_shape = (32, 48))

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            movie_view = window.movie_view
            visible = movie_view.mapToScene(movie_view.viewport().rect()).boundingRect()
            image_rect = QtCore.QRectF(0.0, 0.0, 48.0, 32.0)
            assert (visible.contains(image_rect))

            # Fit means filling the view in one direction, not just being inside it.
            assert ((visible.width() < 1.1*image_rect.width()) or (visible.height() < 1.1*image_rect.height()))
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_capture():
    """
    Frame capture saves the movie view as a PNG.
    """
    movie_name, h5_name = makeTestData("test_visualizer_capture", twoNeighbors())
    capture_name = storm_analysis.getPathOutputTest("test_visualizer_capture.png")
    storm_analysis.removeFile(capture_name)

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            window.loadMovie(movie_name)
            window.captureToFile(capture_name)
            image = QtGui.QImage(capture_name)
            pixel_ratio = window.movie_view.devicePixelRatioF()
            assert (not image.isNull())
            assert (image.width() == round(window.movie_view.width()*pixel_ratio))
            assert (image.height() == round(window.movie_view.height()*pixel_ratio))
            assert (window.reported_errors == [])
        finally:
            window.close()

    assert (exceptions == [])


def test_visualizer_contrast_slider_range():
    """
    The contrast slider runs from 0 to the smallest power of 2, and at
    least 256, that covers the brightest pixel shown so far and the
    maximum setting. It grows as needed and starts again with each movie.
    """
    dim_frame = numpy.full((32, 32), 100.0)
    dim_frame[2, 3] = 900.0
    bright_frame = dim_frame.copy()
    bright_frame[5, 5] = 5000.0
    flat_frame = numpy.full((32, 32), 100.0)

    movie_names = {}
    for name, frames in [["bright", [dim_frame, bright_frame]],
                         ["dim", [dim_frame]],
                         ["flat", [flat_frame]]]:
        movie_names[name] = storm_analysis.getPathOutputTest("test_visualizer_" + name + ".dax")
        movie = datawriter.DaxWriter(movie_names[name])
        for frame in frames:
            movie.addFrame(frame)
        movie.close()

    def sliderTop(window):
        return window.rangeSlider.start + window.rangeSlider.scale

    with RecordExceptions() as exceptions:
        window = makeWindow()
        try:
            # Before any movie the default maximum of 2000 sets the range.
            assert (window.rangeSlider.start == 0)
            assert (sliderTop(window) == 2048)

            # It grows for a brighter frame, and doesn't shrink again.
            window.loadMovie(movie_names["bright"])
            assert (sliderTop(window) == 2048)
            window.incCurFrame(1)
            assert (sliderTop(window) == 8192)
            window.incCurFrame(-1)
            assert (sliderTop(window) == 8192)

            # A maximum above the range extends it.
            window.ui.maxSpinBox.setValue(10000)
            assert (sliderTop(window) == 16384)
            assert (window.rangeSlider.getValues() == [100, 10000])

            # A new movie starts again.
            window.ui.maxSpinBox.setValue(300)
            window.loadMovie(movie_names["dim"])
            assert (sliderTop(window) == 1024)

            window.ui.minSpinBox.setValue(0)
            window.ui.maxSpinBox.setValue(150)
            window.loadMovie(movie_names["flat"])
            assert (sliderTop(window) == 256)
            assert (window.rangeSlider.getValues() == [0, 150])
        finally:
            window.close()

    assert (exceptions == [])


if (__name__ == "__main__"):
    test_visualizer_load_and_step()
    test_visualizer_click_selects_localization()
    test_visualizer_click_after_contrast_change()
    test_visualizer_equal_contrast_limits()
    test_visualizer_bad_files_are_reported()
    test_visualizer_new_movie_closes_old_one()
    test_visualizer_integer_fields()
    test_visualizer_movie_fits_view()
    test_visualizer_capture()
    test_visualizer_contrast_slider_range()
