param([Parameter(Mandatory=$true)][string]$Dataset)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$matlab = 'D:\Matlab R2019a\bin\matlab.exe'
$pilotSeed = 20260717
$values = @(1,10,100)
$scoreRows = @()

foreach ($alpha in $values) {
    foreach ($beta in $values) {
        $tag = ('a{0}_b{1}' -f $alpha,$beta).Replace('.','p')
        $out = Join-Path $root "pilot_snmfwlp_${Dataset}_$tag"
        $log = Join-Path $root "pilot_snmfwlp_${Dataset}_$tag.log"
        $raw = Join-Path $out 'SNMFWLP_common_stop_1seed_raw.csv'
        if (-not (Test-Path -LiteralPath $raw)) {
            $env:SNMFWLP_DATASET = $Dataset
            $env:SNMFWLP_OUTPUT = $out
            $env:SNMFWLP_ALPHA = [string]$alpha
            $env:SNMFWLP_BETA = [string]$beta
            $env:SNMFWLP_SEED_START = [string]$pilotSeed
            $args = @('-nojvm','-nodesktop','-nosplash','-singleCompThread',
                '-logfile',$log,'-sd',$root,'-r','driver_snmfwlp_pilot_env')
            $p = Start-Process -FilePath $matlab -ArgumentList $args -WindowStyle Hidden -PassThru
            $deadline = (Get-Date).AddHours(2)
            while ((Get-Date) -lt $deadline -and -not (Test-Path -LiteralPath $raw)) {
                Start-Sleep -Seconds 10
            }
            if (-not (Test-Path -LiteralPath $raw)) {
                throw "SNMFWLP pilot timed out: dataset=$Dataset alpha=$alpha beta=$beta"
            }
            Get-CimInstance Win32_Process -Filter "Name='MATLAB.exe'" -ErrorAction SilentlyContinue |
                Where-Object { $_.CommandLine -like "*$log*" } |
                ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
        }
        $r = Import-Csv -LiteralPath $raw | Select-Object -First 1
        $scoreRows += [pscustomobject]@{
            dataset=$Dataset; pilot_seed=$pilotSeed; alpha=$alpha; beta=$beta
            validation_ACC=[double]$r.validation_ACC
            iterations=[int]$r.iterations
            stop_residual=[double]$r.stop_residual
            converged=[int]$r.converged
            objective_nonincrease=[int]$r.objective_nonincrease
        }
    }
}

$scorePath = Join-Path $root "pilot_snmfwlp_${Dataset}_scores.csv"
$scoreRows | Sort-Object @{Expression='validation_ACC';Descending=$true},alpha,beta |
    Export-Csv -LiteralPath $scorePath -NoTypeInformation -Encoding UTF8
$best = $scoreRows | Sort-Object @{Expression='validation_ACC';Descending=$true},alpha,beta | Select-Object -First 1
"BEST dataset=$Dataset alpha=$($best.alpha) beta=$($best.beta) validation_ACC=$($best.validation_ACC)" |
    Set-Content -LiteralPath (Join-Path $root "pilot_snmfwlp_${Dataset}_best.txt") -Encoding UTF8
Write-Output $best
