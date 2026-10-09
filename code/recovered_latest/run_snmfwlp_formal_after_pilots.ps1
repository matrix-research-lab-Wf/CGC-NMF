$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$matlab = 'D:\Matlab R2019a\bin\matlab.exe'
$pilotAudit = Join-Path $root 'snmfwlp_pilot_master.log'
$formalAudit = Join-Path $root 'snmfwlp_formal_master.log'

while (-not ((Test-Path -LiteralPath $pilotAudit) -and
        (Select-String -LiteralPath $pilotAudit -Pattern '^ALL_COMPLETE ' -Quiet))) {
    Start-Sleep -Seconds 30
}

foreach ($dataset in @('PIE','YaleB','COIL20','Optdigits','MNIST','COIL100')) {
    $scorePath = Join-Path $root "pilot_snmfwlp_${dataset}_scores.csv"
    if (-not (Test-Path -LiteralPath $scorePath)) { throw "Missing $scorePath" }
    $best = Import-Csv -LiteralPath $scorePath |
        Where-Object { [double]$_.alpha -in @(1,10,100) -and [double]$_.beta -in @(1,10,100) } |
        Sort-Object @{Expression={[double]$_.validation_ACC};Descending=$true},@{Expression={[int]$_.converged};Descending=$true},@{Expression={[int]$_.iterations};Ascending=$true},@{Expression={[double]$_.alpha};Ascending=$true},@{Expression={[double]$_.beta};Ascending=$true} |
        Select-Object -First 1
    $out = Join-Path $root "formal_snmfwlp_3seed_$dataset"
    $log = Join-Path $root "formal_snmfwlp_3seed_$dataset.log"
    $raw = Join-Path $out 'SNMFWLP_common_stop_3seed_raw.csv'
    $decision = Join-Path $out 'SNMFWLP_common_stop_3seed_decision.txt'
    $complete = (Test-Path -LiteralPath $decision) -and (Test-Path -LiteralPath $raw) -and (@(Import-Csv -LiteralPath $raw).Count -eq 3)
    if (-not $complete) {
        $env:SNMFWLP_DATASET = $dataset
        $env:SNMFWLP_OUTPUT = $out
        $env:SNMFWLP_ALPHA = [string]$best.alpha
        $env:SNMFWLP_BETA = [string]$best.beta
        Add-Content -LiteralPath $formalAudit -Value ("START $dataset alpha=$($best.alpha) beta=$($best.beta) $((Get-Date).ToString('o'))")
        $args = @('-nojvm','-nodesktop','-nosplash','-singleCompThread','-logfile',$log,'-sd',$root,'-r','driver_snmfwlp_formal_env')
        $p = Start-Process -FilePath $matlab -ArgumentList $args -WindowStyle Hidden -PassThru
        $deadline = (Get-Date).AddHours(6)
        $complete = $false
        while ((Get-Date) -lt $deadline -and -not $complete) {
            Start-Sleep -Seconds 15
            $complete = (Test-Path -LiteralPath $decision) -and (Test-Path -LiteralPath $raw) -and (@(Import-Csv -LiteralPath $raw).Count -eq 3)
        }
        if (-not $complete) { throw "Formal SNMFWLP timed out for $dataset" }
        Get-CimInstance Win32_Process -Filter "Name='MATLAB.exe'" -ErrorAction SilentlyContinue |
            Where-Object { $_.CommandLine -like "*$log*" } |
            ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
        Add-Content -LiteralPath $formalAudit -Value ("END $dataset $((Get-Date).ToString('o'))")
    }
}
Add-Content -LiteralPath $formalAudit -Value ("ALL_COMPLETE $((Get-Date).ToString('o'))")
