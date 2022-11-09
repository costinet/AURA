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


//get requested action
if (array_key_exists("action",$_GET)) {
	$action = htmlspecialchars($_GET["action"]);
} else {
	die('Error: No action specified');
}


//Handshake request
if (strcmp($action, "handshake") == 0) {
	echo ini_get('post_max_size').', '.ini_get('upload_max_filesize');
	return;
} 

// Uploading data
elseif (strcmp($action, "upload") == 0  || strcmp($action, "sync") == 0) {
	//get and decode data
	$DB = $_POST["data"];
	$DB = json_decode($DB, FALSE);  // FALSE makes it return as an object
	
	if (json_last_error() <> JSON_ERROR_NONE) {
		die("Error: Invalid JSON formatting of data");
	} 
	
	$numFETS = count($DB);

	//For each FET in the user database
	if ($numFETS >1) {
		foreach ($DB as $fet) {
			$paramData = getTransistorParams($fet, $userName);
			//syncParamData($paramData, $conn);
			
			syncParamData($conn, $fet, $userName);
			
			$graphData = getTransistorGraphs($fet, $userName);
			syncGraphData($conn, $graphData, $paramData, $userName);
		}
	} else {
		$fet = $DB;
		$paramData = getTransistorParams($fet, $userName);
		$graphData = getTransistorGraphs($fet, $userName);
		
		//echo var_dump($graphData[0]);
		
		//echo 'Number of existing graphs:'.checkGraphDuplicates($graphData[0], $paramData, $userName, $conn).PHP_EOL;
		//echo var_dump($graphData);
		
		syncParamData($conn, $fet, $userName);
		syncGraphData($conn, $graphData, $paramData, $userName);
		
	}
	
	//updateUserLastSync($conn, $userName);
	
} 

elseif (strcmp($action, "download") == 0 || strcmp($action, "sync") == 0) {
	//echo "download<br>";
	//$lastSyncTime = strtotime(getUserLastSync($conn, $userName)." -1 month");   // "-1 month" ONLY FOR DEBUGGING
	$lastSyncTime = strtotime(getUserLastSync($conn, $userName)); 
	//echo "lastSync Unixtime: ".$lastSyncTime."<br>";
	//echo "lastSync Readable: ".date("Y-m-d",$lastSyncTime)."<br>";
	
	$newGraphs = getGraphDataAfterDate($conn, $lastSyncTime);
	if (!empty($newGraphs)) {
		echo "<graphData>";
		echo $newGraphs;
		echo "<\graphData><br>".PHP_EOL;
	} else {
		echo "No New Graphs<br>";
	}
	
	$newParams= getParamsDataAfterDate($conn, $lastSyncTime);
	if (!empty($newParams)) {
		echo "<paramData>";
		echo $newParams;
		echo "<\paramData><br>".PHP_EOL;
	} else {
		echo "No New Parameters<br>";
	}
	
	//updateUserLastSync($conn, $userName);
	
}

// Invalid action
else {
	die("invalid action specified");
}


$conn->close();
echo "This ends our interation for today.";
	
	