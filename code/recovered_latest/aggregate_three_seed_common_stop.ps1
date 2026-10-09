$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$datasets = @('PIE','YaleB','COIL20','Optdigits','MNIST')
$rows = @()

foreach ($dataset in $datasets) {
    $path = Join-Path $root "gate3_core_$dataset\common_stop_core_3seed_raw.csv"
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing $path" }
    $rows += Import-Csv -LiteralPath $path
}

foreach ($pair in @(
    @('gate3_gnmfld_COIL100','common_stop_gnmfld_3seed_raw.csv'),
    @('gate3_gocpair_COIL100','common_stop_goc_pair_3seed_raw.csv')
)) {
    $path = Join-Path (Join-Path $root $pair[0]) $pair[1]
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing $path" }
    $rows += Import-Csv -LiteralPath $path
}

$methodOrder = @{'GNMFLD'=1;'GOCNMF'=2;'CGC-GOCNMF'=3}
$datasetOrder = @{'PIE'=1;'YaleB'=2;'COIL20'=3;'COIL100'=4;'Optdigits'=5;'MNIST'=6}
$rows = $rows | Sort-Object @{Expression={[int]$datasetOrder[$_.dataset]}},@{Expression={[int]$_.seed}},@{Expression={[int]$methodOrder[$_.method]}}

if ($rows.Count -ne 54) { throw "Expected 54 rows, found $($rows.Count)." }
$keys = $rows | ForEach-Object { "$($_.dataset)|$($_.seed)|$($_.method)" }
if (($keys | Sort-Object -Unique).Count -ne 54) { throw 'Duplicate dataset/seed/method keys.' }

$badConvergence = @($rows | Where-Object { [int]$_.converged -ne 1 })
$badObjective = @($rows | Where-Object { [int]$_.objective_nonincrease -ne 1 })
$badResidual = @($rows | Where-Object { [double]$_.stop_residual -gt 1e-4 })

$rawOut = Join-Path $root 'common_stop_3seed_all_methods_raw.csv'
$summaryOut = Join-Path $root 'common_stop_3seed_all_methods_summary.csv'
$auditOut = Join-Path $root 'common_stop_3seed_all_methods_audit.txt'
$rows | Export-Csv -LiteralPath $rawOut -NoTypeInformation -Encoding UTF8

$summary = foreach ($group in ($rows | Group-Object dataset,method)) {
    $g = @($group.Group)
    $acc = @($g | ForEach-Object {[double]$_.unlabeled_ACC})
    $nmi = @($g | ForEach-Object {[double]$_.unlabeled_NMI})
    $iter = @($g | ForEach-Object {[double]$_.iterations_used})
    $accMean = ($acc | Measure-Object -Average).Average
    $nmiMean = ($nmi | Measure-Object -Average).Average
    $iterMean = ($iter | Measure-Object -Average).Average
    $accSd = [math]::Sqrt((($acc | ForEach-Object {($_-$accMean)*($_-$accMean)} | Measure-Object -Sum).Sum)/($acc.Count-1))
    $nmiSd = [math]::Sqrt((($nmi | ForEach-Object {($_-$nmiMean)*($_-$nmiMean)} | Measure-Object -Sum).Sum)/($nmi.Count-1))
    [pscustomobject]@{dataset=$g[0].dataset;method=$g[0].method;runs=$g.Count;ACC_mean=$accMean;ACC_sd=$accSd;NMI_mean=$nmiMean;NMI_sd=$nmiSd;iterations_mean=$iterMean;converged="$(@($g | Where-Object {[int]$_.converged -eq 1}).Count)/$($g.Count)";objective_nonincrease="$(@($g | Where-Object {[int]$_.objective_nonincrease -eq 1}).Count)/$($g.Count)"}
}
$summary | Sort-Object @{Expression={[int]$datasetOrder[$_.dataset]}},@{Expression={[int]$methodOrder[$_.method]}} | Export-Csv -LiteralPath $summaryOut -NoTypeInformation -Encoding UTF8

$pass = ($badConvergence.Count -eq 0 -and $badObjective.Count -eq 0 -and $badResidual.Count -eq 0)
@(
    'UNIFIED COMMON-STOP 3-SEED AGGREGATE AUDIT'
    'criterion=row-normalized representation relative change <= 1e-4 for 5 consecutive iterations after minIter=20; maxIter=1000'
    "rows=$($rows.Count)/54"
    "convergence_failures=$($badConvergence.Count)"
    "objective_nonincrease_failures=$($badObjective.Count)"
    "terminal_residual_failures=$($badResidual.Count)"
    "UNIFIED_COMMON_STOP_3SEED_PASS=$([int]$pass)"
) | Set-Content -LiteralPath $auditOut -Encoding UTF8

if (-not $pass) { throw "Aggregate audit failed. See $auditOut" }
Write-Output "UNIFIED_COMMON_STOP_3SEED_PASS=1"
Write-Output $rawOut
Write-Output $summaryOut
Write-Output $auditOut
