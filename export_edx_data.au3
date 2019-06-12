#include <GUIConstantsEx.au3>
#include <MsgBoxConstants.au3>

; ===============================================================================================================================
; Title .........: AZtec EDX map exporter
; Description ...: Simple windows utility to export EDX maps from the AZtec EDX software
; Author(s) .....: Daniel Martin Lundeby
; Date ..........: May 28th, 2019
; ===============================================================================================================================

; Create main window
$hMain = GUICreate("Export EDX maps from AZtec", 400, 200)
GUISetState(@SW_SHOW)
Opt("GUIOnEventMode", 1)
GUISetOnEvent($GUI_EVENT_CLOSE, "CloseButton")

; Create folder selection button
$path = ""
$iChooseFolder = GUICtrlCreateButton ( "Select folder", 10, 30)
GUICtrlSetOnEvent($iChooseFolder, "SelectFolderButton")
$iLabel = GUICtrlCreateLabel("<Folder is not chosen>", 80, 30, 300, 30)

; Create checkboxes
$iTranspose = GUICtrlCreateCheckbox ( "Transpose", 10, 60 )
$iFlipx = GUICtrlCreateCheckbox ( "Flip x", 10, 80 )
$iFlipy = GUICtrlCreateCheckbox ( "Flip y", 10, 100 )

; Create export button
$iExport = GUICtrlCreateButton ( "Export", 10, 130)
GUICtrlSetOnEvent($iExport, "ExportButton")

; Create help button
$iHelp = GUICtrlCreateButton( "Help", 100, 130)
GUICtrlSetOnEvent($iHelp, "OpenHelpFile")


While 1
    Sleep(100) ; Sleep to reduce CPU usage
WEnd

Func CloseButton()
    Exit
EndFunc


Func SelectFolderButton()
	$path = FileSelectFolder ( "Choose folder with files exported from AZtec", "")
	$len = StringLen ( $path)

	$showpath = ""
	$maxlen = 50
	If $len > $maxlen Then
		$showpath = "..." & StringMid($path, $len-$maxlen)
	Else
		$showpath = $path
	EndIf
	GUICtrlSetData($iLabel, $showpath)
EndFunc

Func ExportButton()

	$cmd = "C:/ProgramData/Anaconda3/python.exe ./export_edx_data.py "  & """" & $path & """"

	If GUICtrlRead($iTranspose) == $GUI_CHECKED Then
		$cmd = $cmd & " --transpose"
	EndIf
	If GUICtrlRead($iFlipx) == $GUI_CHECKED Then
		$cmd = $cmd & " --flipx"
	EndIf
	If GUICtrlRead($iFlipy) == $GUI_CHECKED Then
		$cmd = $cmd & " --flipy"
	EndIf

	$program = Run($cmd)

EndFunc

Func OpenHelpFile()
	ShellExecute("Exporting EDX data from AZtec.pdf")
EndFunc