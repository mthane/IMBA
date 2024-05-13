import sys
import cv2
import os
import numpy as np
import datetime
import subprocess
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from distutils.filelist import FileList

LRVTRACK = './lrvTrack.AppImage'
#LRVTRACK = './IMBtracker_cpp/lrvTrack'
LRVTRACK = os.path.abspath(LRVTRACK)
TMPDIR = '/tmp'
#from TrackUI import Ui_TrackExperiment
from TrackUI import Ui_MainWindow

class trackWorker(QtCore.QObject):
    # Signal to signify move to new video
    process_next = pyqtSignal(int, int, name="proccessNext")
    process_end = pyqtSignal()

    #def __init__(self):
    #    super().__init__()

    def setup(self, videos_model, bigDishRadioVal,
                 brightnessVal, contrastVal, gammaVal,
                 collisionVal,
                 ROIx, ROIy, ROIr,
                 outputPath):
        self.videos_model = videos_model
        self.bigDishRadioVal = bigDishRadioVal
        self.contrastVal = contrastVal
        self.collisionVal = collisionVal
        self.gammaVal = gammaVal
        self.brightnessVal = brightnessVal
        self.ROIx = ROIx
        self.ROIy = ROIy
        self.ROIr = ROIr
        self.outputPath = outputPath

        print("Inintializing worker with " + str(self.videos_model.rowCount()))

    #@QtCore.pyqtSlot()
    def processVideos(self):
        print("Starting processing...")
        for row in range(self.videos_model.rowCount()):
            self.currentTrackingID = row
            print("Processing video: " + str(row))
            #print("Name: " + str(self.videos_model.item(row,1).text()))
            self.currentTrackingID = row
            ret = self.processSingleVideo(row)
            self.process_next.emit(row,ret)
            #app.processEvents()
        self.process_end.emit()

    def processSingleVideo(self, ctid):
        if self.videos_model.item(ctid,0).checkState() == QtCore.Qt.Checked:
            video = self.videos_model.item(ctid,1).text()
            trial_name = os.path.basename(os.path.dirname(video))
            rec_name = os.path.basename(os.path.dirname(
                                    os.path.dirname(
                                    video)))
            group_name = os.path.basename(os.path.dirname(
                                            os.path.dirname(
                                            os.path.dirname(
                                            video))))
            exp_name = os.path.basename(os.path.dirname(
                                        os.path.dirname(
                                        os.path.dirname(
                                        os.path.dirname(
                                        video)))))
            subfolder = group_name + '/' + rec_name + '/' + trial_name
            video_path = self.videos_model.item(ctid,1).text()
            trial_path = os.path.dirname(video_path)
            dish_size = 84
            min_size = 110
            if self.bigDishRadioVal:
                dish_size = 138
            b = self.brightnessVal
            c = self.contrastVal
            g = self.gammaVal
            if b == 0.0:
                blist = []
            else:
                blist = ['--brightness'] + [str(b+255)]

            if c == 1.0:
                clist = []
            else:
                clist = ['--contrast'] + [str(c)]

            if g == 1.0:
                glist = []
            else:
                glist = ['--gamma'] + [str(g)]

            if self.ROIr == 0:
                ROIlist = []
            else:
                ROIlist = ["--roi"] + [str(self.ROIx) + "," + str(self.ROIy) + ',' + str(self.ROIr)]

            if self.collisionVal:
                collisionlist = ['--use-model', '--model-duration', '8']
            else:
                collisionlist = []

            print("start tracking")
            print(trial_path)
            trial_path = trial_path.replace("\\", "/")
            
            video_path = video_path.replace("\\", "/")
            print(trial_path)
            command = ' '.join([LRVTRACK,
                        '-x', # offline background computation
                        '-z', # dish size set below
                        str(dish_size)] +
                        blist + clist + glist +
                        ['--min-obj-size', #min size of blob
                        str(min_size),
                        '--max-obj-size', #max size of blob
                        str(5000)] +
                        # ROI if specified
                        ROIlist +
                        collisionlist + [
                        '-p', #run multiprocessor
                        '--thread-count', # number of threads...
                        '9',
                        '-o',  # generate output
                        '-v',
                        '16',
                        '-d',
                        '13',
                        '-w',
                        '0.04',
                        '-u',
                        '-t',
                        '--metadata-file',
                        '-i',
                        video_path]
                        )
            print(command)
            with open(trial_path+"/stdout.log","wb") as out, open(trial_path+"/stderr.log","wb") as err:
                process = subprocess.Popen(command, 
                                            shell=True,
                                            cwd = trial_path,
                                            stderr=err,
                                            stdout=out)
                exitCode = process.wait()
            #trackP.waitForFinished(-1)
            #exitCode = trackP.exitCode()
            print('ExitCode: ' + str(exitCode))
            #ctid = self.currentTrackingID
            output_path = self.outputPath
            if exitCode==0:
                print("test")
                print(output_path+"/"+subfolder)
                if not os.path.exists(output_path + '/' + subfolder):
                    os.makedirs(output_path + '/' + subfolder)

                vidsAndLogsPath = output_path + '/' + subfolder
                real_output = output_path + '/' + subfolder
                # Copy files to output folder
                try:
                    #print("Metadata:")
                    #print("")
                    print("copy metadata")
                    cm = subprocess.call(['cp ' + trial_path + '/metadata.txt ' + vidsAndLogsPath],shell=True)
                    print("move stdout.log")
                    co = subprocess.call(['mv ' + trial_path + '/stdout.log ' + vidsAndLogsPath],shell=True)
                    print("move stderr.log")
                    ce = subprocess.call(['mv ' + trial_path + '/stderr.log ' + vidsAndLogsPath],shell=True)
                    cd = subprocess.call(['mv ' + trial_path + '/*-data/* ' + real_output],shell=True)
                    cv = subprocess.call(['mv ' + trial_path + '/2*.mp4 '  + vidsAndLogsPath],shell=True)
                    
                    cc = subprocess.call(['bash createCollisionsCSVinExperiment.sh ' + vidsAndLogsPath + '/stdout.log'],shell=True)
                except:
                    print("Error: Some of the files were not found")
                    return 1
                return 0
            else:
                return 1
                    #self.videos_model.item(ctid,6).setText('Fail')
            #self.experimentInfoGroupBox.setEnabled(True)
            #self.trackingParamsGroupBox.setEnabled(True)
            #self.trackStartBtn.setEnabled(True)
            #self.trackCancelBtn.setEnabled(False)

            #ctid += 1
            #self.trackProgressBar.setValue(
            #    (ctid*100.0)/self.videos_model.rowCount())
            #self.currentTrackingID = ctid
            #if self.videos_model.rowCount() <= ctid:
            #    return
            #self.trackVidID(ctid)


class trackApp(QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
            # It sets up layout and widgets that are defined
        self.folderBrowseBtn.clicked.connect(self.selectExperimentFolder)
        self.btnSelectAll.clicked.connect(self.selectAll)
        self.btnSelectNone.clicked.connect(self.selectNone)
        self.trackStartBtn.clicked.connect(self.startedTracking)
        self.outputFolderBrowseBtn.clicked.connect(self.selectOutputFolder)
        self.brightnessVal.valueChanged.connect(self.imgProcChange)
        self.contrastVal.valueChanged.connect(self.imgProcChange)
        self.gammaVal.valueChanged.connect(self.imgProcChange)
        self.brightnessCheckbox.stateChanged.connect(self.imgProcChange)
        self.contrastCheckbox.stateChanged.connect(self.imgProcChange)
        self.gammaCheckbox.stateChanged.connect(self.imgProcChange)
        self.zoomCheckbox.stateChanged.connect(self.imgProcChange)
        self.ROIcheckbox.stateChanged.connect(self.showFrame)
        self.ROIXspinbox.valueChanged.connect(self.updateFrame)
        self.ROIYspinbox.valueChanged.connect(self.updateFrame)
        self.ROIRspinbox.valueChanged.connect(self.updateFrame)
        self.videos_model = QStandardItemModel(self)
        self.outputAppendCheckBox.setVisible(False)
        self.tableView.setModel(self.videos_model)
        self.tableView.setEditTriggers(QAbstractItemView.NoEditTriggers)
        track_check_header = QStandardItem('Track')
        track_check_header.setCheckable(True)
        track_check_header.setCheckState(2)
        self.videos_model.setHorizontalHeaderItem(0,track_check_header)
        self.videos_model.setHorizontalHeaderItem(1,QStandardItem('path_str'))
        self.videos_model.setHorizontalHeaderItem(2,QStandardItem('Experiment'))
        self.videos_model.setHorizontalHeaderItem(3,QStandardItem('Group'))
        self.videos_model.setHorizontalHeaderItem(4,QStandardItem('Reciprocal'))
        self.videos_model.setHorizontalHeaderItem(5,QStandardItem('Trial'))
        self.videos_model.setHorizontalHeaderItem(6,QStandardItem('Result'))
        #self.tableView.horizontalHeader().setStretchLastSection(True)
        self.tableView.horizontalHeader().setSectionResizeMode(5,QHeaderView.Stretch)
        self.tableView.horizontalHeader().resizeSection(0,40)
        self.tableView.horizontalHeader().resizeSection(6,40)
        self.tableView.setColumnHidden(1,True)
        self.tableView.verticalHeader().setVisible(False)
        self.tableView.clicked.connect(self.fullRowSelect)
        self.tableView.setSelectionBehavior(QTableView.SelectRows)
        self.tableView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.trackProgressBar.setFormat('%v of %m done')
        self.trackProgressBar.setMaximum(1)
        # set current tracking id to -1
        self.currentTrackingID = -1
        self.cv_image = None
        self.v_thread = QtCore.QThread()
        self.videos_track = trackWorker()
        # Settings for specification of ROI
        self.ROI_r = 0
        self.ROI_x = 0
        self.ROI_y = 0
        self.resultsTabWidget.setCurrentIndex(1)
        #time.sleep(1)
        #self.resultsTabWidget.setCurrentIndex(0)

    def updateFrame(self, int):
        if not self.ROIcheckbox.checkState():
            self.ROI_x = 0
            self.ROI_y = 0
            self.ROI_r = 0
        else:
            self.ROI_x = self.ROIXspinbox.value()
            self.ROI_y = self.ROIYspinbox.value()
            self.ROI_r = self.ROIRspinbox.value()
        self.showFrame()
        self.imgProcChange()
        self.showFrame()

    def selectAll(self):
        for row in range(self.videos_model.rowCount()):
            if self.videos_model.item(row,0).checkState() != QtCore.Qt.Checked:
                self.videos_model.item(row,0).setCheckState(QtCore.Qt.Checked)
    def selectNone(self):
        for row in range(self.videos_model.rowCount()):
            if self.videos_model.item(row,0).checkState() != QtCore.Qt.Unchecked:
                self.videos_model.item(row,0).setCheckState(QtCore.Qt.Unchecked)

    def showFrame(self):
        # Our operations on the frame come here
        #image = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        cv_image = self.cv_image
        print(cv_image)
        h, w, _ = cv_image.shape
        image = QImage(cv_image, w, h, 3 * w,
                                 QImage.Format_RGB888)
        if self.ROIRspinbox.value() == 0:
            #print("Setting initial ROI values")
            minDIM = min(cv_image.shape[0],cv_image.shape[1])
            self.ROI_x = cv_image.shape[0]/2
            self.ROI_y = cv_image.shape[1]/2
            self.ROI_r = 0.9 * minDIM/2
            self.ROIXspinbox.setValue(int(self.ROI_x))
            self.ROIYspinbox.setValue(int(self.ROI_y))
            self.ROIRspinbox.setValue(int(self.ROI_r))
        pix = QPixmap(image)
        if self.zoomCheckbox.isChecked():
            self.framePixmap = pix
            if self.ROIcheckbox.checkState():
                ROICircle = QPainter()
                pen = QPen(QtCore.Qt.red, 2, QtCore.Qt.SolidLine)
                ROICircle.begin(pix)
                ROICircle.setPen(pen)
                ROICircle.drawEllipse(QtCore.QPoint(
                                        self.ROIXspinbox.value(),
                                        self.ROIYspinbox.value()),
                                      self.ROIRspinbox.value(),
                                      self.ROIRspinbox.value())
                ROICircle.end()
            self.frameShowLbl.setPixmap(self.framePixmap)
        else:
            w = self.frameShowScrollArea.width()
            h = self.frameShowScrollArea.height()
            #print(self.frameShowScrollArea.sizeHint())
            v = w
            if w>h:
                v = h
            if self.ROIcheckbox.checkState():
                ROICircle = QPainter()
                pen = QPen(QtCore.Qt.red, 4, QtCore.Qt.SolidLine)
                ROICircle.begin(pix)
                ROICircle.setPen(pen)
                ROICircle.drawEllipse(QtCore.QPoint(
                                        self.ROIXspinbox.value(),
                                        self.ROIYspinbox.value()),
                                      self.ROIRspinbox.value(),
                                      self.ROIRspinbox.value())
                ROICircle.end()
            npix = pix.scaled(v,v,QtCore.Qt.KeepAspectRatio)
            self.framePixmap = npix
            self.frameShowLbl.setPixmap(self.framePixmap)

    def fullRowSelect(self,selection):
        if len(self.tableView.selectionModel().selectedIndexes()) > 0:
            vid_path = self.videos_model.item(selection.row(),1).text()
            cap = cv2.VideoCapture(vid_path)
            ret, cv_image= cap.read()
            self.cv_image = cv_image
            self.cv_frame = cv_image.copy()
            self.showFrame()
        else:
            self.cv_image = None
            self.frameShowLbl.clear()
            self.frameShowLbl.repaint()

    def selectOutputFolder(self):
        path = QFileDialog.getExistingDirectory(
            options=QFileDialog.ShowDirsOnly)
        self.outputPathLineEdit.setText(path)

    def selectExperimentFolder(self):
        path = QFileDialog.getExistingDirectory(
            options=QFileDialog.ShowDirsOnly)
        self.pathLineEdit.setText(path)
        if path == "":
            return
        try:
            self.outputPathLineEdit.setText(path + "/Tracked_" + str(datetime.date.today()))
            exp_search = FileList()
            exp_search.findall(path)
            exp_search.include_pattern("*.mp4", anchor=0)
            self.experiment_videos = exp_search.files
            self.loadVideos()
        except:
            path = ""
            self.outputPathLineEdit.setText("")
            self.experiment_videos = []
            self.load_videos()

    def imgProc(self,b,c,g):
        o = self.cv_frame.copy()
        if c!=1.0 or b!=0.0:
            #b = b-255
            o = cv2.addWeighted(o,0,o,c,b,o)
        if g!=1.0:
            o = np.float32(o)
            #o = o/255.0
            np.clip(o, 0.0, 1.0, out=o)
            cv2.pow(o,g,o)
            o = o * 255.0
            o = np.uint8(o)
            np.clip(o, 0, 255, out=o)
        meanVal = np.mean(o)
        o[o<meanVal] = meanVal
        cv2.normalize(o,o,0,255,cv2.NORM_MINMAX)
        self.cv_image = o
        self.showFrame()

    def imgProcChange(self):
        if self.cv_frame is None:
            self.frameShowLbl.clear()
            self.frameShowLbl.repaint()
        else:
            c = 1.0
            b = 0.0
            g = 1.0
            if self.brightnessCheckbox.isChecked():
                b = self.brightnessVal.value()
            if self.contrastCheckbox.isChecked():
                c = self.contrastVal.value()
            if self.gammaCheckbox.isChecked():
                g = self.gammaVal.value()
            self.imgProc(b,c,g)

    def loadVideos(self):
        #Clear treeview before filling
        #Also clear image processing stuff before as well
        if self.videos_model.hasChildren():
            self.videos_model.removeRows(0,self.videos_model.rowCount())
            self.cv_image = None
            self.frameShowLbl.clear()
            self.frameShowLbl.repaint()
            self.imgProcChange()

        for video in self.experiment_videos:
            trial_name = QStandardItem(
            os.path.basename(os.path.dirname(video)))
            print(trial_name)
            print(video)
            rec_name = QStandardItem(
                os.path.basename(os.path.dirname(
                                 os.path.dirname(
                                 video))))
            group_name = QStandardItem(
                os.path.basename(os.path.dirname(
                                          os.path.dirname(
                                          os.path.dirname(
                                          video)))))
            exp_name = QStandardItem(
                os.path.basename(os.path.dirname(
                                        os.path.dirname(
                                        os.path.dirname(
                                        os.path.dirname(
                                        video))))))
            track_check = QStandardItem('')
            track_check.setCheckable(True)
            track_check.setCheckState(2)
            self.videos_model.appendRow([track_check,
                                         QStandardItem(video),
                                         exp_name,
                                         group_name,
                                         rec_name,
                                         trial_name,
                                         QStandardItem('')
                                         ])
        self.trackProgressBar.setMaximum(self.videos_model.rowCount())

    def onNewVid(self,idx):
        total = self.videos_model.rowCount()
        self.trackProgressBar.setValue(idx)

    def endedTracking(self):
        self.trackProgressBar.setValue(self.trackProgressBar.maximum())
        self.experimentInfoGroupBox.setEnabled(True)
        self.trackingParamsGroupBox.setEnabled(True)
        self.trackStartBtn.setEnabled(True)
        self.trackCancelBtn.setEnabled(False)
        self.v_thread.quit()
        self.v_thread.wait()

    def updateStatus(self, item_no, result):
        if result == 0:
            val = 'OK'
        else:
            val = "Fail"
        self.videos_model.item(item_no,6).setText(val)
        self.onNewVid(item_no)

    def startedTracking(self):
        self.trackProgressBar.setValue(0)
        self.experimentInfoGroupBox.setEnabled(False)
        self.trackingParamsGroupBox.setEnabled(False)
        self.trackStartBtn.setEnabled(False)
        self.trackCancelBtn.setEnabled(True)
        if self.bigDishRadioBtn.isChecked():
            self.bigDishRadioVal = True
        else:
            self.bigDishRadioVal = False
        if self.withColResRadioBtn.isChecked():
            self.collisionVal = True
        else:
            self.collisionVal = False

        if self.brightnessCheckbox.isChecked():
            self.brightness = self.brightnessVal.value()
        else:
            self.brightness = 0.0

        if self.contrastCheckbox.isChecked():
            self.contrast = self.contrastVal.value()
        else:
            self.contrast = 1.0

        if self.gammaCheckbox.isChecked():
            self.gamma = self.gammaVal.value()
        else:
            self.gamma = 1.0

        print("Set up trackworker")
        print("ROI X:" + str(self.ROI_x))
        self.videos_track = trackWorker()
        self.videos_track.setup(self.videos_model,
                               self.bigDishRadioVal,
                               self.brightness,
                               self.contrast,
                               self.gamma,
                               self.collisionVal,
                               self.ROI_x,
                               self.ROI_y,
                               self.ROI_r,
                               self.outputPathLineEdit.text())
        self.videos_track.moveToThread(self.v_thread)
        print("moved to thread")
        self.videos_track.process_end.connect(self.v_thread.quit)
        self.videos_track.process_next.connect(self.updateStatus)
        self.v_thread.started.connect(self.videos_track.processVideos)
        self.v_thread.finished.connect(self.endedTracking)
        self.v_thread.start()
        self.videos_track.deleteLater()
        print("Started thread")


def main():
    app = QApplication(sys.argv)  # A new instance of QApplication
    form = trackApp()                 # We set the form to be our ExampleApp (design)
    form.show()                         # Show the form
    app.exec_()                         # and execute the app


if __name__ == '__main__':              # if we're running file directly and not importing it
    main()                              # run the main function
