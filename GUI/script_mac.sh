#!/bin/sh

mkdir BayeScEnv_GUI/Contents/Frameworks

cp -R /opt/local/Library/Frameworks/QtCore.framework/ BayeScEnv_GUI.app/Contents/Frameworks/QtCore.framework

cp -R /opt/local/Library/Frameworks/QtGui.framework/ BayeScEnv_GUI.app/Contents/Frameworks/QtGui.framework

cp -R /opt/local/lib/libz.1.2.8.dylib BayeScEnv_GUI.app/Contents/Frameworks/libz.1.dylib

cp -R /opt/local/lib/libpng16.16.dylib BayeScEnv_GUI.app/Contents/Frameworks/libpng16.16.dylib

install_name_tool -change /opt/local/Library/Frameworks/QtGui.framework/Versions/4/QtGui @executable_path/../Frameworks/QtGui.framework/Versions/4/QtGui BayeScEnv_GUI.app/Contents/MacOS/BayeScEnv_GUI

install_name_tool -change /opt/local/Library/Frameworks/QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore.framework/Versions/4/QtCore BayeScEnv_GUI.app/Contents/MacOS/BayeScEnv_GUI

install_name_tool -change /opt/local/Library/Frameworks/QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore.framework/Versions/4/QtCore BayeScEnv_GUI.app/Contents/Frameworks/QtGui.framework/Versions/4/QtGui

install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/../Frameworks/libz.1.2.8.dylib BayeScEnv_GUI.app/Contents/Frameworks/QtCore.framework/Versions/4/QtCore

install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/../Frameworks/libz.1.2.8.dylib BayeScEnv_GUI.app/Contents/Frameworks/QtGui.framework/Versions/4/QtGui

install_name_tool -change /opt/local/lib/libpng16.16.dylib @executable_path/../Frameworks/libpng16.16.dylib BayeScEnv_GUI.app/Contents/Frameworks/QtGui.framework/Versions/4/QtGui