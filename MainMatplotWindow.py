# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\MainMatplotWindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.

import matplotlib
import matplotlib.pyplot as pyplot

matplotlib.use('Qt5Agg')

from PyQt5 import QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from MainWindow import Ui_MainWindow
from Modules import *

import sys, os, os.path

class MainMatplotWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainMatplotWindow, self).__init__()

        self.userInterface = Ui_MainWindow()
        self.userInterface.setupUi(self)

        self.figure = pyplot.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.userInterface.verticalLayout.addWidget(self.toolbar)
        self.userInterface.verticalLayout.addWidget(self.canvas)

        icondir = ""
        if hasattr(sys, "_MEIPASS"): icondir = os.path.join(sys._MEIPASS, './icons/mainIcon.png')
        else: icondir = './icons/mainIcon.png'
        self.setWindowIcon(QtGui.QIcon(icondir))

        self.setHandlers()
        self.setValidators()

        self.dataContext = DataContext()
        self.resultSolution = None
    
    def setHandlers(self):
        self.userInterface.ClearSheduleButton.clicked.connect(self.onClickClearSheduleButton)
        self.userInterface.RunCalculateButton.clicked.connect(self.onClickRunCalculateButton)
        self.userInterface.ShowPhiFunctionButton.clicked.connect(self.onClickShowPhiFunctionButton)

    def setValidators(self):
        self.doubleValidator = QtGui.QDoubleValidator()
        self.userInterface.StringLengthLineEdit.setValidator(self.doubleValidator)
        self.userInterface.TimeLineEdit.setValidator(self.doubleValidator)
        self.userInterface.LengthStepLineEdit.setValidator(self.doubleValidator)
        self.userInterface.TimeStepLineEdit.setValidator(self.doubleValidator)
        self.userInterface.CoefficientALineEdit.setValidator(self.doubleValidator)
        self.userInterface.FZeroLineEdit.setValidator(self.doubleValidator)
        self.userInterface.FFirstLineEdit.setValidator(self.doubleValidator)
        self.userInterface.FSecondLineEdit.setValidator(self.doubleValidator)
        self.userInterface.FThirdLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PhiZeroLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PhiFirstLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PhiSecondLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PsiZeroLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PsiFirtsLineEdit.setValidator(self.doubleValidator)
        self.userInterface.PsiSecondLineEdit.setValidator(self.doubleValidator)
    
    def onClickClearSheduleButton(self):
        self.resultSolution = None
        pyplot.gcf().clear()
        self.canvas.draw()

    def onClickRunCalculateButton(self):
        if(not self.fillDataContext()): return

        k = (self.dataContext.coefficientA * self.dataContext.coefficientA * self.dataContext.timeStep) / (self.dataContext.lengthStep * self.dataContext.lengthStep)
        if (k < 0.25):
            self.resultSolution = Solution(self.dataContext)
            self.drawShedule(self.resultSolution)
        else: ShowErorMassageBox("Шаги по длине и по времени не удовлетворяют условию устойчивости:\na*a*τ/(h*h) < 1/4")

    def onClickShowPhiFunctionButton(self):
        if(self.resultSolution == None):
            ShowErorMassageBox("Вывести функцию Φ нельзя вывести, так как холст пуст!")
            return
        
        pyplot.plot(self.resultSolution[0], self.resultSolution[3], 'g')
        self.canvas.draw()
    
    def drawShedule(self, data):
        pyplot.gcf().clear()
        self.canvas.draw()

        pyplot.grid()
        pyplot.xlabel("Длина струны")
        pyplot.ylabel("Время")
        pyplot.plot(data[0], data[1], 'b')

        self.canvas.draw() 

    def fillDataContext(self):
        if self.userInterface.StringLengthLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.TimeLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.LengthStepLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.TimeStepLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.CoefficientALineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.FZeroLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.FFirstLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.FSecondLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.FThirdLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PhiZeroLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PhiFirstLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PhiSecondLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PsiZeroLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PsiFirtsLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        if self.userInterface.PsiSecondLineEdit.text() == "":
            ShowErorMassageBox("Заполните все поля для ввода!")
            return False
        
        self.dataContext.initialize(
            float(self.userInterface.StringLengthLineEdit.text().replace(',','.')),
            float(self.userInterface.TimeLineEdit.text().replace(',','.')),
            float(self.userInterface.LengthStepLineEdit.text().replace(',','.')),
            float(self.userInterface.TimeStepLineEdit.text().replace(',','.')),
            float(self.userInterface.CoefficientALineEdit.text().replace(',','.')),
            float(self.userInterface.FZeroLineEdit.text().replace(',','.')),
            float(self.userInterface.FFirstLineEdit.text().replace(',','.')),
            float(self.userInterface.FSecondLineEdit.text().replace(',','.')),
            float(self.userInterface.FThirdLineEdit.text().replace(',','.')),
            float(self.userInterface.PhiZeroLineEdit.text().replace(',','.')),
            float(self.userInterface.PhiFirstLineEdit.text().replace(',','.')),
            float(self.userInterface.PhiSecondLineEdit.text().replace(',','.')),
            float(self.userInterface.PsiZeroLineEdit.text().replace(',','.')),
            float(self.userInterface.PsiFirtsLineEdit.text().replace(',','.')),
            float(self.userInterface.PsiSecondLineEdit.text().replace(',','.'))
        )

        return True