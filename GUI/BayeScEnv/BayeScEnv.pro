#-------------------------------------------------
#
# Project created by QtCreator 2014-08-21T13:04:29
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BayeScEnv_GUI
TEMPLATE = app


SOURCES += main.cpp\
        gui.cpp

HEADERS  += gui.h

FORMS    += gui.ui

CONFIG+=/openmp

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp -lpthread -static-libgcc -static-libstdc++
