from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMessageBox

import math
import numpy
import sys, os, os.path

def ShowErorMassageBox(text):
    msg = QMessageBox()
    msg.setWindowTitle("Ошибка")
    icondir = ""
    if hasattr(sys, "_MEIPASS"): icondir = os.path.join(sys._MEIPASS, './icons/warning.png')
    else: icondir = './icons/warning.png'
    msg.setWindowIcon(QIcon(icondir))
    msg.setInformativeText(text)
    msg.setStyleSheet("min-width: 200px;")
    msg.exec_()

def ShowInformationMassageBox(text):
    msg = QMessageBox()
    msg.setWindowTitle("Информация")
    icondir = ""
    if hasattr(sys, "_MEIPASS"): icondir = os.path.join(sys._MEIPASS, './icons/info.png')
    else: icondir = './icons/info.png'
    msg.setWindowIcon(QIcon(icondir))
    msg.setInformativeText(text)
    msg.setStyleSheet("min-width: 200px;")
    msg.exec_()

# Class for convenient data storage
class DataContext:
    def __init__(self):
        self.initialized = False

        self.stringLength = None
        self.time = None
        self.lengthStep = None
        self.timeStep = None

        self.coefficientA = None

        self.f0 = None
        self.f1 = None
        self.f2 = None
        self.f3 = None

        self.phi0 = None
        self.phi1 = None
        self.phi2 = None
        self.phi3 = None

        self.psi0 = None
        self.psi1 = None
        self.psi2 = None
        self.psi3 = None

        self.MainFunction = None
        self.PhiFunction = None
        self.PsiFunction = None
    
    def initialize(self, stringLength, time, lengthStep, timeStep, coefficientA, f0, f1, f2, f3, phi0, phi1, phi2, psi0, psi1, psi2):
        self.stringLength = stringLength
        self.time = time
        self.lengthStep = lengthStep
        self.timeStep = timeStep

        self.coefficientA = coefficientA

        self.f0 = f0
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3

        self.phi0 = phi0
        self.phi1 = phi1
        self.phi2 = phi2

        self.psi0 = psi0
        self.psi1 = psi1
        self.psi2 = psi2

        self.phi3 = 1 - (stringLength * phi0 + stringLength * 0.5 * phi1 + stringLength * 0.5 * phi2)
        self.psi3 = (0 - (stringLength * phi0 * psi0 + stringLength * 0.5 * phi1 * psi1 + stringLength * 0.5 * phi2 * psi2)) / self.phi3

        self.MainFunction = lambda x, t: self.f0 + self.f1 * math.cos((math.pi * x) / self.stringLength) + self.f2 * math.cos((2 * math.pi * x) / self.stringLength) + self.f3 * math.cos((3 * math.pi * x) / self.stringLength)
        self.PhiFunction = lambda x: self.phi0 + self.phi1 * math.cos((math.pi * x) / self.stringLength) + self.phi2 * math.cos((2 * math.pi * x) / self.stringLength) + self.phi3 * math.cos((3 * math.pi * x) / self.stringLength)
        self.PsiFunction = lambda x: self.psi0 + self.psi1 * math.cos((math.pi * x) / self.stringLength) + self.psi2 * math.cos((2 * math.pi * x) / self.stringLength) + self.psi3 * math.cos((3 * math.pi * x) / self.stringLength)

        self.initialized = True

# Solve u_tt = a^2 * u_xx + f(x,t) on (0 < x < l) and (0 <= t <= T)
def Solution(data):
    count_T = int(data.time / data.timeStep)
    t = numpy.linspace(0, count_T * data.timeStep, count_T + 1)

    count_N = int(data.stringLength / (data.timeStep * data.coefficientA / float(data.lengthStep)))
    x = numpy.linspace(0, data.stringLength, count_N + 1)

    dx = x[1] - x[0]
    dt = t[1] - t[0]

    phiFunction = numpy.zeros(count_N + 1)
    u     = numpy.zeros(count_N + 1)
    u_n   = numpy.zeros(count_N + 1)
    u_nm1 = numpy.zeros(count_N + 1)

    for i in range(0, count_N + 1):
        phiFunction[i] = data.PhiFunction(x[i])
        u_n[i] = data.PhiFunction(x[i])

    n = 0
    for i in range(1, count_N):
        u[i] = u_n[i] + dt * data.PsiFunction(x[i]) + 0.5 * data.lengthStep**2 * (u_n[i-1] - 2 * u_n[i] + u_n[i+1]) + 0.5 * dt**2 * data.MainFunction(x[i], t[n])

    u[0] = 0
    u[count_N] = 0

    u_nm1[:] = u_n
    u_n[:] = u

    for n in range(1, count_T):
        for i in range(1, count_N):
            u[i] = - u_nm1[i] + 2 * u_n[i] + data.lengthStep**2 * (u_n[i-1] - 2 * u_n[i] + u_n[i+1]) + dt**2 * data.MainFunction(x[i], t[n])

        u[0] = 0
        u[count_N] = 0

        u_nm1[:] = u_n
        u_n[:] = u

    return x, u, t, phiFunction
