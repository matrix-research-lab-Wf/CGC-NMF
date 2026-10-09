$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$resultRoot = Join-Path $root 'results\recovered_six_dataset_3seed'
$datasets = @('PIE','YaleB','COIL20','Optdigits','MNIST','COIL100')

$common = @()
$snm = @()
$audit = @('RECOVERED SIX-DATASET THREE-SEED AUDIT')
foreach ($dataset in $datasets) {
    $commonPath = Join-Path $resultRoot "common_stop_$dataset\common_stop_core_3seed_summary.csv"
    $snmPath = Join-Path $resultRoot "snmfwlp_$dataset\SNMFWLP_common_stop_3seed_summary.csv"
    $common += Import-Csv -LiteralPath $commonPath
    $snm += Import-Csv -LiteralPath $snmPath

    $commonDecision = Join-Path $resultRoot "common_stop_$dataset\common_stop_core_3seed_decision.txt"
    $snmDecision = Join-Path $resultRoot "snmfwlp_$dataset\SNMFWLP_common_stop_3seed_decision.txt"
    $commonPass = Select-String -LiteralPath $commonDecision -Pattern '^THREE_SEED_INTEGRITY_PASS=1$' -Quiet
    $snmFinite = Select-String -LiteralPath $snmDecision -Pattern '^FINITE_OUTPUT_PASS=1$' -Quiet
    $raw = Import-Csv (Join-Path $resultRoot "snmfwlp_$dataset\SNMFWLP_common_stop_3seed_raw.csv")
    $snmConverged = @($raw | Where-Object { $_.converged -eq '1' }).Count
    $snmObjective = @($raw | Where-Object { $_.objective_nonincrease -eq '1' }).Count
    $audit += ("{0}: common_integrity={1}; snm_finite={2}; snm_converged={3}/3; snm_objective_nonincrease={4}/3" -f `
        $dataset,[int]$commonPass,[int]$snmFinite,$snmConverged,$snmObjective)
}

$common | Export-Csv -LiteralPath (Join-Path $root 'results\recovered_common_stop_3seed_summary.csv') -NoTypeInformation -Encoding UTF8
$snm | Export-Csv -LiteralPath (Join-Path $root 'results\recovered_snmfwlp_3seed_summary.csv') -NoTypeInformation -Encoding UTF8
$audit | Set-Content -LiteralPath (Join-Path $root 'results\recovered_3seed_audit.txt') -Encoding UTF8

$gitRoot = Join-Path $root '.git'
$manifest = Get-ChildItem -LiteralPath $root -Recurse -File -Force |
    Where-Object {
        $_.Name -ne 'MANIFEST_SHA256.csv' -and
        -not $_.FullName.StartsWith($gitRoot, [System.StringComparison]::OrdinalIgnoreCase)
    } |
    ForEach-Object {
        $relative = $_.FullName.Substring($root.Length + 1)
        [pscustomobject]@{
            path = $relative
            bytes = $_.Length
            sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        }
    }
$manifest | Sort-Object path | Export-Csv -LiteralPath (Join-Path $root 'MANIFEST_SHA256.csv') -NoTypeInformation -Encoding UTF8
