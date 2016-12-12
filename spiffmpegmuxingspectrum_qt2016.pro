TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += spiffmpegmuxingspectrum_qt2016.c

INCLUDEPATH += "..\lib-src\ffmpeg\ffmpeg-20160113-git-9ca64c3-win32-dev\ffmpeg-20160113-git-9ca64c3-win32-dev\include"
INCLUDEPATH += "..\lib-src\libsndfile\include"
INCLUDEPATH += "..\lib-src\rfftw_qt2016"
INCLUDEPATH += "..\lib-src\freeimage\Source"

DEFINES += __STDC_CONSTANT_MACROS
DEFINES += __STDC_FORMAT_MACROS

LIBS += -L"..\lib-src\ffmpeg\ffmpeg-20160113-git-9ca64c3-win32-dev\ffmpeg-20160113-git-9ca64c3-win32-dev\lib" -lavcodec -lavformat -lswscale -lavutil -lswresample
LIBS += -L"..\lib-src\libsndfile(rel)" -llibsndfile-1
LIBS += -L"..\lib-src\rfftw_qt2016-release\release" -lrfftw
LIBS += -L"..\lib-src\freeimage\Release" -lFreeImage
