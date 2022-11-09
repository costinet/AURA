<?php

ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

require('E:\Dropbox\UTK\Research\AURA\AURAdb\@netConnect\webServer\matlabInterface.php');

$tableParamNames = getTransistorTableColNames();
$conn = powerLibConnect();

//validate user API key
if (array_key_exists("apiKey",$_GET)) {
	$key = htmlspecialchars($_GET["apiKey"]);
	$key = mysqli_real_escape_string($conn, $key);
} else {
	die('Error: No API key specified');
}
$userName = validateUser($conn, $key);


//echo "download<br>";
//$lastSyncTime = strtotime(getUserLastSync($conn, $userName)." -1 month");   // "-1 month" ONLY FOR DEBUGGING
//$lastSyncTime = strtotime(getUserLastSync($conn, $userName)); 
//echo "lastSync Unixtime: ".$lastSyncTime."<br>";
//echo "lastSync Readable: ".date("Y-m-d",$lastSyncTime)."<br>";

$lastSyncTime = strtotime( "2000-01-31");
$newGraphs = getGraphDataAfterDate($conn, $lastSyncTime);
$newParams= getParamsDataAfterDate($conn, $lastSyncTime);


$paramArray = json_decode($newParams);
$graphArray = json_decode($newGraphs);

$uniqeDevices = array();
foreach ($paramArray as $param){
	$uniqeDevices = array_merge($uniqeDevices, array($param->partNumber));
}
foreach ($graphArray as $param){
	$uniqeDevices = array_merge($uniqeDevices, array($param->partNumber));
}
$uniqeDevices = array_unique($uniqeDevices);

// Display Library

$page = file_get_contents('libView.php');
$page = str_replace('##Styles##',file_get_contents('styles.html'),$page);
$page = str_replace('##Scripts##',file_get_contents('scripts.html'),$page);


$FETtemplate = file_get_contents('templates/FETformat.txt');
$FETtable = '';
foreach ($uniqeDevices as $device){
	$newEntry = str_replace('##PartNo##',$device,$FETtemplate);
	
	$devParams = '';
	foreach ($paramArray as $param){
		if (strcmp($param->partNumber,$device) == 0) {
			$devParams  = var_export($param, true);
		}
	}
	
	$devGraphs = '';
	foreach ($graphArray as $graph){
		if (strcmp($graph->partNumber,$device) == 0) {
			$devGraphs  = $devGraphs.$graph->title.'<br>';
		}
	}
	
	$newEntry = str_replace('##Params##','Params:<br>'.$devParams.'<br>Graphs:<br>'.$devGraphs,$newEntry);
	
	$FETtable = $FETtable.PHP_EOL.$newEntry;
}

$page = str_replace('##FET_TABLE##',$FETtable,$page);

echo($page);

//echo var_dump($graphArray);





$conn->close();
echo "This ends our interation for now.";
	
	