# Lanza scripts/calib_sigma.py DESACOPLADO de la sesión SSH (Windows).
# En Windows, lo que se lanza dentro de una sesión de OpenSSH muere al desconectarse;
# un proceso creado por WMI (Win32_Process.Create) no, así que sigue corriendo.
#
#   powershell -ExecutionPolicy Bypass -File scripts\launch_detached.ps1                 # corrida completa
#   powershell -ExecutionPolicy Bypass -File scripts\launch_detached.ps1 --workers 12    # con argumentos
#
# Ver progreso:  Get-Content results\calib_sigma\log.txt -Tail 5 -Wait -Encoding UTF8
# Cortar:        powershell -ExecutionPolicy Bypass -File scripts\stop_detached.ps1
# Sin bloque param: todo lo que se pase va tal cual a calib_sigma.py.
$ExtraArgs = $args -join " "

$root = Split-Path -Parent $PSScriptRoot
$py = Join-Path $root ".venv\Scripts\python.exe"
$outDir = Join-Path $root "results\calib_sigma"
$log = Join-Path $outDir "log.txt"
New-Item -ItemType Directory -Force $outDir | Out-Null

$cmd = "cmd /c `"`"$py`" -u scripts\calib_sigma.py $ExtraArgs >> `"$log`" 2>&1`""
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create `
    -Arguments @{ CommandLine = $cmd; CurrentDirectory = $root }
if ($r.ReturnValue -ne 0) { throw "No se pudo lanzar (ReturnValue=$($r.ReturnValue))" }
Set-Content -Path (Join-Path $outDir "pid.txt") -Value $r.ProcessId
"Lanzado (PID $($r.ProcessId)). Log: $log"
"Progreso: Get-Content `"$log`" -Tail 5 -Wait -Encoding UTF8"
"Argumentos: $ExtraArgs"
