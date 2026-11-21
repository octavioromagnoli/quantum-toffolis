# Corta la corrida lanzada con launch_detached.ps1 (mata el árbol de procesos).
# Lo ya calculado queda en results\calib_sigma\runs.jsonl; relanzar retoma desde ahí.
$root = Split-Path -Parent $PSScriptRoot
$pidFile = Join-Path $root "results\calib_sigma\pid.txt"
if (-not (Test-Path $pidFile)) { "No hay pid.txt"; exit }
$p = (Get-Content $pidFile).Trim()
taskkill /PID $p /T /F
Remove-Item $pidFile
