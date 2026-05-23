$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$RawDir = Join-Path $ProjectRoot "data\raw"
$ExtractDir = Join-Path $RawDir "GSE165722_RAW"

New-Item -ItemType Directory -Force -Path $RawDir | Out-Null
New-Item -ItemType Directory -Force -Path $ExtractDir | Out-Null

$RawTar = Join-Path $RawDir "GSE165722_RAW.tar"
$SeriesMatrix = Join-Path $RawDir "GSE165722_series_matrix.txt.gz"

$RawTarUrl = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165722/suppl/GSE165722_RAW.tar"
$SeriesMatrixUrl = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165722/matrix/GSE165722_series_matrix.txt.gz"

if (-not (Test-Path $RawTar)) {
    Write-Host "Downloading GSE165722_RAW.tar..."
    curl.exe --ssl-no-revoke -L $RawTarUrl -o $RawTar
} else {
    Write-Host "GSE165722_RAW.tar already exists. Skipping download."
}

if (-not (Test-Path $SeriesMatrix)) {
    Write-Host "Downloading GSE165722_series_matrix.txt.gz..."
    curl.exe --ssl-no-revoke -L $SeriesMatrixUrl -o $SeriesMatrix
} else {
    Write-Host "GSE165722_series_matrix.txt.gz already exists. Skipping download."
}

Write-Host "Extracting raw tar..."
tar.exe -xf $RawTar -C $ExtractDir

Write-Host "Done."
Write-Host "Raw files are in: $ExtractDir"
