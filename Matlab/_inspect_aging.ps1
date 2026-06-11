$path = 'd:\NARP_Data\RTrack_NARPMale\TgF344_Aging_AllRats.xlsx'
$excel = New-Object -ComObject Excel.Application
$excel.Visible = $false
$wb = $excel.Workbooks.Open($path)
foreach ($ws in $wb.Sheets) {
    Write-Output "--- $($ws.Name) ---"
    $headers = @()
    for ($c = 1; $c -le 20; $c++) {
        $v = $ws.Cells.Item(1, $c).Text
        if ($v) { $headers += $v }
    }
    Write-Output ("Columns: " + ($headers -join ', '))
    for ($r = 2; $r -le 6; $r++) {
        $row = @()
        for ($c = 1; $c -le [Math]::Min(8, $headers.Count); $c++) {
            $row += $ws.Cells.Item($r, $c).Text
        }
        Write-Output ($row -join ' | ')
    }
}
$wb.Close($false)
$excel.Quit()
[System.Runtime.Interopservices.Marshal]::ReleaseComObject($excel) | Out-Null
