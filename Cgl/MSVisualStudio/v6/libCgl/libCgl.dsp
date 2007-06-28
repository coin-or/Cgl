# Microsoft Developer Studio Project File - Name="libCgl" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libCgl - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libCgl.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libCgl.mak" CFG="libCgl - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libCgl - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libCgl - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libCgl - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GR /GX /O2 /I "..\..\..\..\Cgl\src\CglDuplicateRow" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding2" /I "..\..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\..\Cgl\src\CglClique" /I "..\..\..\..\Cgl\src\CglOddHole" /I "..\..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\..\Cgl\src\CglGomory" /I "..\..\..\..\Cgl\src\CglPreProcess" /I "..\..\..\..\Cgl\src\CglProbing" /I "..\..\..\..\Cgl\src" /I "..\..\..\..\Osi\src" /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libCgl - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "..\..\..\..\Cgl\src\CglDuplicateRow" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding2" /I "..\..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\..\Cgl\src\CglClique" /I "..\..\..\..\Cgl\src\CglOddHole" /I "..\..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\..\Cgl\src\CglGomory" /I "..\..\..\..\Cgl\src\CglPreProcess" /I "..\..\..\..\Cgl\src\CglProbing" /I "..\..\..\..\Cgl\src" /I "..\..\..\..\Osi\src" /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libCgl - Win32 Release"
# Name "libCgl - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglAllDifferent\CglAllDifferent.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglClique\CglClique.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglClique\CglCliqueHelper.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglCutGenerator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglDuplicateRow\CglDuplicateRow.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglFlowCover\CglFlowCover.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglGomory\CglGomory.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglKnapsackCover\CglKnapsackCover.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglLiftAndProject\CglLiftAndProject.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMessage.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding\CglMixedIntegerRounding.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding2\CglMixedIntegerRounding2.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglOddHole\CglOddHole.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglPreProcess\CglPreProcess.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglProbing\CglProbing.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglRedSplit\CglRedSplit.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglSimpleRounding\CglSimpleRounding.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglStored.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglTwomir\CglTwomir.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglAllDifferent\CglAllDifferent.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglClique\CglClique.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglConfig.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglCutGenerator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglDuplicateRow\CglDuplicateRow.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglFlowCover\CglFlowCover.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglGomory\CglGomory.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglKnapsackCover\CglKnapsackCover.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglLiftAndProject\CglLiftAndProject.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding\CglMixedIntegerRounding.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding2\CglMixedIntegerRounding2.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglOddHole\CglOddHole.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglPreProcess\CglPreProcess.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglProbing\CglProbing.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglRedSplit\CglRedSplit.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglSimpleRounding\CglSimpleRounding.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglStored.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglTwomir\CglTwomir.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinBuild.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinDistance.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinError.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFactorization.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFinite.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFloatEqual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinHelperFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinIndexedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessageHandler.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVectorBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPragma.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinShallowPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinSort.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinTime.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStart.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStartBasis.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiColCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCollections.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCuts.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiPresolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiRowCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiRowCutDebugger.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiSolverInterface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiSolverParameters.hpp
# End Source File
# End Group
# End Target
# End Project
