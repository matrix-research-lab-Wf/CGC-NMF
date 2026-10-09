$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$matlab = 'D:\Matlab R2019a\bin\matlab.exe'
$resultRoot = Join-Path $root 'results'
$jobs = @(
    @{ Dataset='YaleB'; Alpha='10'; Beta='100' },
    @{ Dataset='COIL20'; Alpha='10'; Beta='100' },
    @{ Dataset='Optdigits'; Alpha='10'; Beta='100' },
    @{ Dataset='MNIST'; Alpha='10'; Beta='100' },
    @{ Dataset='COIL100'; Alpha='10'; Beta='100' }
)
foreach ($job in $jobs) {
    $env:RECOVERED_DATASET = $job.Dataset
    $env:RECOVERED_SNM_ALPHA = $job.Alpha
    $env:RECOVERED_SNM_BETA = $job.Beta
    $log = Join-Path $resultRoot ("recovered_3seed_{0}.log" -f $job.Dataset)
    $arguments = @('-nodesktop','-nosplash','-singleCompThread',
        '-logfile',$log,'-sd',$PSScriptRoot,'-r','driver_recovered_dataset_env')
    $process = Start-Process -FilePath $matlab -ArgumentList $arguments `
        -WindowStyle Hidden -PassThru -Wait
    if ($process.ExitCode -ne 0) {
        throw "MATLAB failed for $($job.Dataset); inspect $log"
    }
}
