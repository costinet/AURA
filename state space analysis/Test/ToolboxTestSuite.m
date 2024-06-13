classdef ToolboxTestSuite < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        % function unimplementedTest(testCase)
        %     testCase.verifyFail("Unimplemented test");
        % end

        function checkPLECS(testCase)
            plecs()
        end

        function checkComponentAddParam(testCase)
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
        end

        function checkComponentSubsref(testCase)
            t=transistor('testTransistor');
            t.addParameter('Ron', 5 , 'typ');
            t.addParameter('Coss', 5 , 'min');
            t.addParameter('Qg', 5 , 'max'); 
            
            tDB = transistorDB();
            tDB.add(t);

            
            len = length(tDB.components);
            
            assert(isa(tDB.componentType, 'transistor'), 'Failed 0')
            
            testCase.verifyEqual(length(tDB(1:3)),3)
            testCase.verifyEqual(length({tDB(1:3).partNumber}),3)
            testCase.verifyEqual(length(tDB(1).parameters) , length(tDB.components(1).parameters));
            testCase.verifyEqual(length({tDB.components.partNumber}) , len);
            testCase.verifyEqual(length(tDB.graphs), len);
            
            testCase.verifyEqual(length({tDB(1:2).ron.typ}) , 2);
            testCase.verifyEqual(length(tDB(1).ron.typ) , 1);
            testCase.verifyEqual(length({tDB(:).ron.max}) , len);
            
            fet1 = tDB(1);
            param1 = fet1.parameters(1);
            name1 = {param1.name};
            
            testCase.verifyEqual(fet1.parameters(1),param1);
            testCase.verifyEqual(fet1.parameters(1).name , param1(1).name);
            testCase.verifyTrue(strcmp(tDB(1).parameters(1).name(1), 'R'));
            testCase.verifyEqual({fet1.parameters.name} , name1);
            testCase.verifyEqual({tDB(1).parameters(1).name} , name1);
        end

        function checkNetlisParserMethods(testCase)
            %% From SMPSim call
            sim = SMPSim('3LevelBuck.net');
            sim.u = [100 1 1 1 1]';
            Xss1 = sim.steadyState;

            swvec = evalin('base','swvec');
            ts = evalin('base','ts');
            u = [100 1 1 1 1]';


            %% From Netlist File
            sim2 = SMPSim();
            sim2.converter.topology.circuitParser = NetlistCircuitParser(sim2.converter.topology);
            % sim2.topology.loadCircuit('3LevelBuck.net')
            sim2.parser.loadModel('3LevelBuck.net');
            sim2.u = u;
            sim2.converter.setSwitchingPattern(swvec,ts);
            Xss2 = sim2.steadyState;
            
            %% From Test Netlist
            fid=fopen('3LevelBuck.net');
            netlist = char(fread(fid)'); 
            fclose(fid);
            
            sim3 = SMPSim();
            sim3.converter.topology.circuitParser = NetlistCircuitParser(sim3.converter.topology);
            sim3.parser.loadModel(netlist)
            sim3.u = u;
            sim3.converter.setSwitchingPattern(swvec,ts);
            Xss3 = sim3.steadyState;
            
            %% From Stored Topology
            storedTopology = sim3.parser.saveTopology('Three Level Buck');

            sim4 = SMPSim();
            sim4.converter.topology.circuitParser = NetlistCircuitParser(sim4.converter.topology);
            sim4.parser.loadTopology(storedTopology)
            sim4.u = u;
            sim4.converter.setSwitchingPattern(swvec,ts);
            Xss4 = sim4.steadyState;

            testCase.verifyEqual(Xss1,Xss2,1e-6);
            testCase.verifyEqual(Xss1,Xss3,1e-6);
            testCase.verifyEqual(Xss1,Xss4,1e-6);
        end

    end

end