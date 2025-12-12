# 最大同时运行的 MATLAB 数量，根据你电脑情况调整
$maxParallel = 5

# MATLAB 可执行文件路径：请按你自己的版本修改
$matlab = "D:\matlab\bin\matlab.exe"

# 想要扫描的 d13C_CO2 和 d13C_methane 取值
$d13C_CO2_vals     = @(-5, -6, -7, -8, -9)
$d13C_methane_vals = @(-36, -38, -39, -40)

foreach ($co2 in $d13C_CO2_vals) {
    foreach ($ch4 in $d13C_methane_vals) {

        # 如果当前运行的 MATLAB 个数 >= 最大并行数，就等一等
        while ( (Get-Process matlab -ErrorAction SilentlyContinue | Measure-Object).Count -ge $maxParallel ) {

            $count = (Get-Process matlab -ErrorAction SilentlyContinue | Measure-Object).Count
            Write-Host ("[{0}] {1} MATLAB processes running, waiting..." -f (Get-Date -Format HH:mm:ss), $count) -ForegroundColor Yellow

            Start-Sleep -Seconds 10
        }

        # 有空位了，启动一个新任务
        Write-Host ("[{0}] Launching CHEES_circulation({1}, {2})" -f (Get-Date -Format HH:mm:ss), $co2, $ch4) -ForegroundColor Cyan

        Start-Process $matlab `
            -ArgumentList @(
                "-nosplash",
                "-nodesktop",
                "-r",
	"try; cd('$workdir'); CHEES_circulation($co2,$ch4); catch ME; disp(getReport(ME)); end; exit"
             ) `
            -WindowStyle Hidden
    }
}

Write-Host ("[{0}] All jobs submitted." -f (Get-Date -Format HH:mm:ss)) -ForegroundColor Green
