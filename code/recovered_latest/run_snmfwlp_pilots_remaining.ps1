$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$audit = Join-Path $root 'snmfwlp_pilot_master.log'
foreach ($dataset in @('YaleB','COIL20','Optdigits','MNIST','COIL100')) {
    Add-Content -LiteralPath $audit -Value ("START $dataset $((Get-Date).ToString('o'))")
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File (Join-Path $root 'run_snmfwlp_labeled_pilot.ps1') -Dataset $dataset
    if ($LASTEXITCODE -ne 0) { throw "Pilot failed for $dataset" }
    Add-Content -LiteralPath $audit -Value ("END $dataset $((Get-Date).ToString('o'))")
}
Add-Content -LiteralPath $audit -Value ("ALL_COMPLETE $((Get-Date).ToString('o'))")
