%% component->addParameter
if(1)
    % Test 1: valid parameters
    t=transistor;
    % t.addParameter('Inductance', 5)
    try 
        t.addParameter('Inductance', 5)
    catch e 
        err = 1;
    end
    assert(err == 1, 'Failed Test 1');
    
    % Test 2: valid parameters
    t=transistor;
    t.addParameter('Ron', 5);
    
    assert( strcmp(t.parameters(1).name, 'Rds') && ...
        t.parameters(1).typ == 5000 && ...
        all(strcmp(t.parameters(1).unit, {'m'  'Ohm'})), 'Failed Test 2')
    
    % Test 3: triple param
    t=transistor;
    t.addParameter('Ron', [5 10 2]);
    
    assert( strcmp(t.parameters(1).name, 'Rds') && ...
        t.parameters(1).typ == 5000 && ...
        t.parameters(1).max == 10000 && ...
        t.parameters(1).min == 2000 && ...
        all(strcmp(t.parameters(1).unit, {'m'  'Ohm'})), 'Failed Test 3')
    
    % Test 4: test conditions
    t=transistor;
    t.addParameter('Ron', [5 10 2], 'test = 15 V');
    
    assert( strcmp(t.parameters(1).conditions, 'test = 15 V') , 'Failed Test 4')
    
    % Test 5: direct set max min
    t=transistor;
    t.addParameter('Ron', 5 , 'typ');
    t.addParameter('Coss', 5 , 'min');
    t.addParameter('Qg', 5 , 'max');
    
    assert(t.parameters(1).typ == 5000 && ...
        t.parameters(2).min == 5e12 && ...
        t.parameters(3).max == 5e9 , 'Failed Test 5')
    
    % Test 6: EPC 8002
    tDB = transistorDB;
    fet = tDB(1);
    fet.addParameter('Vgs', [0 6 -4]);
    fet.addParameter('Ids', 2, 'max', {'Ta' '=' '25' '°C'});
    fet.addParameter('Idpulse', 2, 'max', {'Ta' '=' '25' '°C'})
    fet.addParameter('Idss', 100e-6, 'max', {'Vds' '=' '52' 'V'})
    
    disp('component->addParameter tests completed successfully');
end

%% component->subsref
if(1)
   t=transistor('testTransistor');
   t.addParameter('Ron', 5 , 'typ');
   t.addParameter('Coss', 5 , 'min');
   t.addParameter('Qg', 5 , 'max'); 
   
   tDB = transistorDB;
   tDB.add(t);
   
   len = length(tDB.components);
   
   assert(isa(tDB.componentType, 'transistor'), 'Failed 0')
   
   assert(length(tDB(1:3)) == 3, 'Failed 1');
   assert(length({tDB(1:3).partNumber}) == 3, 'Failed 2');
   assert(length(tDB(1).parameters) == length(tDB.components(1).parameters), 'Failed 3');
   assert(length({tDB.components.partNumber}) == len, 'Failed 4');
   assert(length({tDB(:).graphs}) == len, 'Failed 5');
   
   assert(length({tDB(1:2).ron.typ}) == 2, 'Failed 6');
   assert(length(tDB(1).ron.typ) == 1, 'Failed 7');
   assert(length({tDB(:).ron.max}) == len, 'Failed 8');
   
   fet1 = tDB(1);
   param1 = fet1.parameters;
   name1 = {param1.name};
   
   assert(isequal(fet1.parameters,param1), 'Failed 9');
   assert(isequal(fet1.parameters(1).name , param1(1).name), 'Failed 10');
   assert(strcmp(tDB(1).parameters(1).name(1), 'R'), 'Failed 11');
   assert(isequal({fet1.parameters.name} , name1), 'Failed 12');
   assert(isequal({tDB(1).parameters.name} , name1), 'Failed 13');
   assert(isequal({tDB(1).parameters(:).name} , name1), 'Failed 14');
   
   
   disp('component->subsref tests completed successfully');

%% componentDB->subsref

    tDB(end-1).graphs(end).plot
   
end

