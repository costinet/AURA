<?php

function getTransistorTableColNames() { 
	//return array("Vds", "Vds_pulse", "Vgs", "Vth", "Rds", "Id", "Id_pulse", "Idss", "Igss", "Vsd", "Is", "Pd", "Tj", "Tstg", "Tc", "Qg", "Qgs", "Qgd", "Qg_th", "Qoss", "Qrr", "Rg", "Yfs", "Eoss", "Eon", "Eoff", "Ciss", "Coss", "Crss", "tr", "tf", "td_off", "td_on", "Rtheta_jc", "Rtheta_jb", "Rtheta_ja");
	
	return array('Vds','Vsd','Vgs','Vth','Ids','Isd','Ig','Idpulse','Ispulse','Igss','Idss','Rds','Rg','Qg','Qgs','Qgd','Qoss','Qrr','Trr','Irr','Coss','Ciss','Crss','Cgd','Cds','Cgs','Eon','Eoff','Eoss','Rthjc','Rthjb','Rthja','Pd','Gfs','Tr','Tdon','Tf','Tdoff','Tj','Tc');
}

function powerLibConnect() {
	$servername = "localhost:3306";
	$username = "root";
	$password = "";
	$dbname = "powerlib";

	$conn = new mysqli($servername, $username, $password, $dbname);
	if ($conn->connect_error) {
		die("Unable to connect to database " . $conn->connect_error);
	} 
	
	//$conn->query("SET NAMES utf8"); 
	//$conn->query("SET CHARACTER SET utf8"); 
	
	return $conn;
}

function updateUserLastSync($conn, $name) {
	
	//not using prepared statement -- $name must come from DB, already, via validateUser
	$query = "UPDATE users SET lastSync=FROM_UNIXTIME(".time().") WHERE name='$name'";
	$result = $conn->query($query);
}

function getUserLastSync($conn, $name) {
	
	//not using prepared statement -- $name must come from DB, already, via validateUser
	$query = "SELECT lastSync FROM users WHERE name='$name'";
	$result = $conn->query($query);
	if ($result->num_rows != 1) {
		echo "Error getting user Last Sync";
		return 0;
	} else {
		$row = $result->fetch_assoc();
		return $row['lastSync'];
	}
}

function validateUser($conn, $key) {
	//$query = "SELECT name FROM users WHERE apiKey='$key'";
	//$result = $conn->query($query);
	
	$stmt = $conn->prepare("SELECT name FROM users WHERE apiKey = ?");
	$stmt->bind_param("s", $key);
	$stmt->execute();
	$result = $stmt->get_result();
	$stmt->close();

	if ($result->num_rows > 0) {
		$row = $result->fetch_assoc();
		return $row['name'];
	} else {
		$conn->close();
		die('Error: specified API key is invalid');
		return 0;
	}
}

function getTransistorParams($fet, $userName) {
	$sqlData = array();
	
	$sqlData['partNumber'] = convertSingleParameter($fet,'partNumber');
	$sqlData['type'] = convertSingleParameter($fet,'type');
	$sqlData['material'] = convertSingleParameter($fet,'material');
	$date = time();
	$sqlData['dateSubmitted'] = "FROM_UNIXTIME($date)";
	$sqlData['submittedBy'] = $userName;
	$sqlData['dataSource'] = 'datasheet';  //THIS NEEDS TO BE CHANGED
	
	// Old (Tyler's) encoding
	$tableParams = convertSingleParameter($fet,'tableParameters');
	if (!empty($tableParams)){
		foreach (getTransistorTableColNames() as $name) {
			if (array_key_exists($name, $tableParams)){
					$param = $tableParams->$name;
					if (array_key_exists('Min', $tableParams->$name)){
						$sqlData[$name.'_min'] = $param->{'Min'};
					} 
					if (array_key_exists('Typ', $tableParams->$name)){
						$sqlData[$name.'_typ'] = $param->{'Typ'};
					} 
					if (array_key_exists('Max', $tableParams->$name)){
						$sqlData[$name.'_max'] = $param->{'Max'};
					} 
				} 
		}
	} else {
	// New (AURA) encoding
		$tableParams = convertSingleParameter($fet,'parameters');
		//foreach (getTransistorTableColNames() as $name) {
		$validParams = getTransistorTableColNames();
		if (!empty($tableParams)){
			foreach ($tableParams as $param) {
				$name =  $param->{'name'};
				if (in_array($name, $validParams)) {
					if (is_numeric($param->{'max'})) {
						$sqlData[$name.'_max'] = $param->{'max'};
					}
					if (is_numeric($param->{'typ'})) {
						$sqlData[$name.'_typ'] = $param->{'typ'};
					}
					if (is_numeric($param->{'min'})) {
						$sqlData[$name.'_min'] = $param->{'min'};
					}
				}
			}
		}
	}
	
	return $sqlData;
}

function getTransistorGraphs($fet, $username) {
	$graphs = convertSingleParameter($fet,'graphs');
	$graphData = array();
	$sqlData = array();
	
	if (!empty($graphs)){
		foreach ($graphs as $graph) {
			$graphData['plotData'] = json_encode($graph->plotData);
			$graphData['testConditions'] = json_encode($graph->testConditions);
			$graphData['title'] = $graph->title;
			$graphData['xLabel'] = $graph->xLabel;
			$graphData['yLabel'] = $graph->yLabel;
			$graphData['dataLabels'] = json_encode($graph->dataLabels);
			
			array_push($sqlData, $graphData);
		}
	}
	
	return $sqlData;
}

function convertSingleParameter($fet, $field) {
	if (is_array($fet)) {
		if (array_key_exists($field, $fet)){
			$val = $fet->{$field};
			if (empty($val)){
				return NULL;
			} else {
				return $val;
			}
		} else {
			return NULL;
		}
	} else {			
		if (property_exists($fet, $field)){
			$val = $fet->{$field};
			if (empty($val)){
				return NULL;
			} else {
				return $val;
			}
		} else {
			return NULL;
		}
	}
	
}

/*function syncParamData($paramData, $conn) {
	// This REPLACE ... COALESCE method only updates values which are NULL, it will not overwrite parameters other than:
	// dateSubmitted, submittedBy,
	$querystring = "REPLACE INTO transistor_table (".PHP_EOL;
	$valuestring = ") VALUES (".PHP_EOL;
	foreach ($paramData as $key => $value) {
		
		$safeKey = mysqli_real_escape_string($conn, $key);
		
		$querystring = $querystring."\t".$safeKey.','.PHP_EOL;
		if (is_numeric($value)) {
			$valuestring = $valuestring."\t COALESCE(".$safeKey.", ".$value.'),'.PHP_EOL;
		} elseif (strcmp($safeKey, 'dateSubmitted')==0) {
			$valuestring = $valuestring."\t".$value.",".PHP_EOL;
		} else {
			//$valuestring = $valuestring."\t'".$value."',".PHP_EOL;
			$valuestring = $valuestring."\t COALESCE('".$value."', ".$safeKey.'),'.PHP_EOL;
		}
	}
	
	$query = substr($querystring, 0, -3).PHP_EOL.substr($valuestring, 0, -3).PHP_EOL.')'.PHP_EOL;

	$result = $conn->query($query);
	if ($result) {
		echo $result."<br>".PHP_EOL;
	} else {
		echo " failed<br>".PHP_EOL;
	}
	
}*/

function syncParamData($conn, $fet, $userName) {
	// This REPLACE ... COALESCE method only updates values which are NULL, it will not overwrite parameters other than:
	// dateSubmitted, submittedBy,
	
	/*$stmt = $conn->prepare("SELECT entryNo FROM transistor_graph WHERE submittedBy = ? AND partNumber = ? AND xLabel = ? AND yLabel = ?");
	$stmt->bind_param("ssss", $userName, $partNumber, $xLabel, $yLabel);
	$stmt->execute();
	$result = $stmt->get_result();
	$stmt->close();*/
	
	//echo var_dump($fet);
	
	$paramData = getTransistorParams($fet, $userName);
	
	$fieldStr = implode(",", array_keys($paramData));
	$vals = array_values($paramData);
	array_splice($vals,array_search("dateSubmitted",array_keys($paramData)),1);
	$typeStr = '';
	$valStr = array();
	foreach ($paramData as $key => $value) {
	//All keys already validated in getTransistorParams
		if (is_numeric($value)) {
			$typeStr = $typeStr."d";
		} elseif (strcmp($key, 'dateSubmitted')==0) {

		} else{
			$typeStr = $typeStr."s";
		}	
		
		if (strcmp($key, 'dateSubmitted')==0) {
			$valStr[] = "NOW()";
		} elseif (strcmp($key, 'partNumber')==0){
			$valStr[] = "?";
		} else {
			$valStr[] = "COALESCE($key, ?)";
		}
	}
	
	$query = "REPLACE INTO transistor_table (".$fieldStr.") VALUES (";
	$query = $query.implode(",",$valStr).")";
	
	
	
	/*echo "New for prepared statement<br>";
	echo $query;
	echo "<br>";
	echo $typeStr;
	echo "<br>";
	echo var_dump(array_merge(array($typeStr), array_values($vals)));
	echo "<br>";*/
	
	$stmt = $conn->prepare($query);
	//$stmt->bind_param($typeStr, ...$vals);
	if (version_compare(phpversion(), '5.6', '<')) {
		$paramArgs = array_merge(array($typeStr), array_values($vals));
		
		$tmp = array();
        foreach($paramArgs as $key => $value) $tmp[$key] = &$paramArgs[$key];
        call_user_func_array(array($stmt, 'bind_param'), $tmp);
		
		//call_user_func_array(array(&$stmt, 'bind_param'), $paramArgs); 
		
		//$ref    = new ReflectionClass('mysqli_stmt');
		//$method = $ref->getMethod("bind_param");
		//$method->invokeArgs($stmt,$paramArgs);
	} else {
		//$stmt->bind_param($typeStr, ...$vals);
		echo("Uncomment the above line if your php version is above 5.6");
		die('PHP version > 5.6 "..." operator needs to be tested');
	}
	//echo var_dump($stmt);
	$stmt->execute();
	//echo var_dump($stmt);
	$result = $stmt->get_result();
	
	
	/*$querystring = "REPLACE INTO transistor_table (".PHP_EOL;
	$valuestring = ") VALUES (".PHP_EOL;
	foreach ($paramData as $key => $value) {
		
		$safeKey = mysqli_real_escape_string($conn, $key);
		
		$querystring = $querystring."\t".$safeKey.','.PHP_EOL;
		if (is_numeric($value)) {
			$valuestring = $valuestring."\t COALESCE(".$safeKey.", ".$value.'),'.PHP_EOL;
		} elseif (strcmp($safeKey, 'dateSubmitted')==0) {
			$valuestring = $valuestring."\t".$value.",".PHP_EOL;
		} else {
			//$valuestring = $valuestring."\t'".$value."',".PHP_EOL;
			$valuestring = $valuestring."\t COALESCE('".$value."', ".$safeKey.'),'.PHP_EOL;
		}
	}
	
	$query = substr($querystring, 0, -3).PHP_EOL.substr($valuestring, 0, -3).PHP_EOL.')'.PHP_EOL;

	$result = $conn->query($query);*/
	
	if ($result) {
		echo $result."<br>".PHP_EOL;
	} else {
		echo $result."<br>";
		echo mysqli_stmt_errno($stmt)."<br>";
		echo " failed<br>".PHP_EOL;
	}
	
	$stmt->close();
	
}

function checkGraphDuplicates($graphData, $paramData, $userName, $conn) {
	
	$partNumber = mysqli_real_escape_string($conn, $paramData['partNumber']);
	$xLabel = mysqli_real_escape_string($conn, $graphData['xLabel']);
	$yLabel = mysqli_real_escape_string($conn, $graphData['yLabel']);
	
	//$querystring = "SELECT entryNo FROM transistor_graph WHERE submittedBy = '$userName' AND partNumber = '$partNumber' AND xLabel='$xLabel' AND yLabel='$yLabel'";
	//$result = $conn->query($querystring);

	$stmt = $conn->prepare("SELECT entryNo FROM transistor_graph WHERE submittedBy = ? AND partNumber = ? AND xLabel = ? AND yLabel = ?");
	$stmt->bind_param("ssss", $userName, $partNumber, $xLabel, $yLabel);
	$stmt->execute();
	$result = $stmt->get_result();
	$stmt->close();
	
	if ($result) {
		return $result->num_rows;
	} else {
		return 0;
	}
}

function syncGraphData($conn, $graphData, $paramData, $userName) {
	$partNumber = mysqli_real_escape_string($conn, $paramData['partNumber']);
	$date = time();
	$dateTime = "FROM_UNIXTIME($date)";
		
	foreach ($graphData as $graph){
		$duplicateData = checkGraphDuplicates($graph, $paramData, $userName, $conn);
		
		/*$xLabel = mysqli_real_escape_string($conn, $graph['xLabel']);
		$yLabel = mysqli_real_escape_string($conn, $graph['yLabel']);
		$title = mysqli_real_escape_string($conn, $graph['title']);
		$testConditions = mysqli_real_escape_string($conn, $graph['testConditions']);
		$dataLabels = mysqli_real_escape_string($conn, $graph['dataLabels']);
		$plotData = mysqli_real_escape_string($conn, $graph['plotData']);*/
		
		$xLabel =  $graph['xLabel'];
		$yLabel =  $graph['yLabel'];
		$title =  $graph['title'];
		$testConditions =  $graph['testConditions'];
		$dataLabels =  $graph['dataLabels'];
		$plotData =  $graph['plotData'];
		
		if ($duplicateData == 0){
			//No duplicates, go ahead and insert
			//echo "No duplicates found for ".$paramData['partNumber']." graph: ".$graph['title'].PHP_EOL;
			
			//$query = "INSERT into transistor_graph (partNumber, dateSubmitted, submittedBy, plotData, testConditions, title, xLabel, yLabel, dataLabels) VALUES (";
			//$query = $query."'$partNumber', $dateTime, '$userName', '$plotData', '$testConditions', '$title', '$xLabel', '$yLabel', '$dataLabels')";
			//$result = $conn->query($query);
			
			/*echo $partNumber.'<br>'.PHP_EOL;
			echo $userName.'<br>'.PHP_EOL;
			echo $plotData.'<br>'.PHP_EOL;
			echo $testConditions.'<br>'.PHP_EOL;
			echo $title.'<br>'.PHP_EOL;
			echo $xLabel.'<br>'.PHP_EOL;
			echo $yLabel.'<br>'.PHP_EOL;
			echo $dataLabels.'<br>'.PHP_EOL;*/
			
			
			$stmt = $conn->prepare("INSERT INTO transistor_graph (partNumber, dateSubmitted, submittedBy, plotData, testConditions, title, xLabel, yLabel, dataLabels) VALUES ( ?, NOW(), ?, ?, ?, ?, ?, ?, ?)");
			//var_dump($stmt);
			$stmt->bind_param("ssssssss", $partNumber, $userName, $plotData, $testConditions, $title, $xLabel, $yLabel, $dataLabels);
			//var_dump($stmt);
			$stmt->execute();
			$result = $stmt->get_result();
			$stmt->close();
			
			if ($result) {
				//echo $query.PHP_EOL;
			} else {
				echo " failed<br>".PHP_EOL;
			}
			
		} else {
			//Duplicates -- maybe insert but flag the old one?  Could still be the same data
			echo "Duplicate found<br>".PHP_EOL;
			echo "I don't know what to do with this, yet";
		}
	}
}

function getGraphDataAfterDate($conn, $lastSyncDate) {
	$query = "SELECT * FROM transistor_graph WHERE dateSubmitted > FROM_UNIXTIME($lastSyncDate)";
	$result = $conn->query($query);

	if ($result->num_rows > 0) {
		$all = $result->fetch_all(MYSQLI_ASSOC);
		return json_encode($all, JSON_UNESCAPED_SLASHES);
	} else {
		return '';
	}
}
	
function getParamsDataAfterDate($conn, $lastSyncDate) {
	$query = "SELECT * FROM transistor_table WHERE dateSubmitted > FROM_UNIXTIME($lastSyncDate)";
	$result = $conn->query($query);

	if ($result->num_rows > 0) {
		$all = $result->fetch_all(MYSQLI_ASSOC);
		if ($result->num_rows > 1) {	
			$i=0;
			foreach ($all as $fet) {
				$fet = array_filter($fet, "notNullFilter");
				$all[$i] = $fet;
				$i = $i+1;
			}
		} else {
			$all = array_filter($all, "notNullFilter");
		}
		
		
		return json_encode($all, JSON_UNESCAPED_SLASHES);
	} else {
		return '';
	}
}

// array_filter function to remove null entries
function notNullFilter($var){
    return ($var !== NULL && $var !== FALSE && $var !== "");
}


	

?>