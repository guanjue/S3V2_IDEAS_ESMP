#-------------------------------------------------
#
# Project created by QtCreator 2017-10-11T02:30:22
#
#-------------------------------------------------

QT       += core gui concurrent

greaterThan(QT_MAJOR_VERSION, 4):

QT += widgets

TARGET = IDEAS_APP
TEMPLATE = app


SOURCES += main.cpp\
        ideas.cpp \
    call_IDEAS.cpp \
    genomicTensor.cpp \
    MixGauss.cpp \
    tensorHMMbase.cpp \
    displaythread.cpp

HEADERS  += ideas.h \
    datastructure.h \
    genomicTensor.h \
    MixGauss.h \
    tensorHMMbase.h \
    displaythread.h

FORMS    += ideas.ui

DISTFILES += \
    tmp.data.tar.gz \
    IDEAS_APP.pro.user \
    makefilestatic \
    t.bed \
    t.input \
    makefile \
    README \
    readme_prepMat \
    readme_runideaspipe

QMAKE_CXXFLAGS += -fopenmp
#QMAKE_CXXFLAGS += -static

LIBS += -fopenmp
#LIBS += -static

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/release/ -lgsl
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/debug/ -lgsl
else:unix: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lgsl

INCLUDEPATH += $$PWD/../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../usr/local/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/release/libgsl.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/debug/libgsl.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/release/gsl.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/debug/gsl.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/libgsl.a

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/release/ -lgsl
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/debug/ -lgsl
else:unix: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lgsl

INCLUDEPATH += $$PWD/../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../usr/local/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/release/ -lgslcblas
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/debug/ -lgslcblas
else:unix: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lgslcblas

INCLUDEPATH += $$PWD/../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../usr/local/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/release/libgslcblas.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/debug/libgslcblas.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/release/gslcblas.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/debug/gslcblas.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/libgslcblas.a

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/release/ -lgslcblas
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/debug/ -lgslcblas
else:unix: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lgslcblas

INCLUDEPATH += $$PWD/../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../usr/local/include
